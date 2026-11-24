#pragma once
#include <numeric>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/ParameterReader.hpp>
#include <GridKit/Model/PhasorDynamics/ControlLimits.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

#include "Exst1Pti.hpp"

namespace GridKit::PhasorDynamics::Exciter
{

  template <class ScalarT, class IdxT>
  Exst1Pti<ScalarT, IdxT>::Exst1Pti(BusT* bus, const DataT& data)
    : monitor_(std::make_unique<MonitorT>(data)), bus_(bus)
  {
    this->size_ = COUNT;
    monitor_->set(DataT::MonitorableVariables::efd, [this]
                  { return this->y_.getData()[EFD]; });
    using Parameters = typename DataT::Parameters;
    const GridKit::Model::ParameterReader reader(data, "Exst1Pti");
    Tr_    = reader.readReal(Parameters::Tr).value_or(Tr_);
    Vimax_ = reader.readReal(Parameters::Vimax).value_or(Vimax_);
    Vimin_ = reader.readReal(Parameters::Vimin).value_or(Vimin_);
    Tc_    = reader.readReal(Parameters::Tc).value_or(Tc_);
    Tb_    = reader.readReal(Parameters::Tb).value_or(Tb_);
    Ka_    = reader.readReal(Parameters::Ka).value_or(Ka_);
    Ta_    = reader.readReal(Parameters::Ta).value_or(Ta_);
    Vrmax_ = reader.readReal(Parameters::Vrmax).value_or(Vrmax_);
    Vrmin_ = reader.readReal(Parameters::Vrmin).value_or(Vrmin_);
    Kc_    = reader.readReal(Parameters::Kc).value_or(Kc_);
    Kf_    = reader.readReal(Parameters::Kf).value_or(Kf_);
    Tf_    = reader.readReal(Parameters::Tf).value_or(Tf_);

    washout_gain_ = Tf_ > 0 ? Kf_ / Tf_ : 0;
    lead_gain_    = Tb_ > 0 ? Tc_ / Tb_ : 0;
  }

  template <class ScalarT, class IdxT>
  Exst1Pti<ScalarT, IdxT>::~Exst1Pti() = default;

  template <class ScalarT, class IdxT>
  int Exst1Pti<ScalarT, IdxT>::allocate()
  {
    if (!this->allocated_)
      this->allocateVectors(COUNT);
    this->tag_.assign(COUNT, false);
    this->variable_indices_.resize(COUNT);
    this->residual_indices_.resize(COUNT);
    std::iota(this->variable_indices_.begin(), this->variable_indices_.end(), IdxT{0});
    std::iota(this->residual_indices_.begin(), this->residual_indices_.end(), IdxT{0});
    this->ws_.resize(Utilities::enum_size<Input>());
    this->ws_indices_.assign(Utilities::enum_size<Input>(), INVALID_INDEX<IdxT>);
    this->wb_.resize(2);
    if (auto& port = ports_.out[Output::efd]; port)
      port.link(this->y_.getData() + EFD, &this->getVariableIndex(EFD));
    this->allocated_ = true;
    return 0;
  }

  template <class ScalarT, class IdxT>
  GridKit::Model::ConfigurationChecks Exst1Pti<ScalarT, IdxT>::verify() const
  {
    GridKit::Model::ConfigurationChecks checks("Exst1Pti");
    checks.check(bus_ != nullptr, "terminal bus required");
    for (auto input : Utilities::enum_values<Input>())
    {
      if (ports_.in[input] && !ports_.in[input].linked())
        checks.fail(std::format("unlinked input {}", Utilities::enum_name(input)));
    }

    checks.check(Ka_ > 0, "Ka must be positive");
    checks.check(Tr_ >= 0 && Ta_ >= 0, "Tr and Ta must be nonnegative");
    checks.check(Tc_ >= 0 && Tb_ >= 0 && (Tb_ > 0 || Tc_ == 0), "Tc requires a positive Tb");
    checks.check(Tf_ >= 0 && (Tf_ > 0 || Kf_ == 0), "Kf requires a positive Tf");
    checks.check(Vimin_ <= Vimax_ && Vrmin_ <= Vrmax_, "voltage limits must be ordered");
    checks.check(Kc_ == 0 || this->ports_.in[Input::ifd], "field current required with commutation loading");
    checks.check(this->ports_.out[Output::efd].connected(), "efd output required");
    return checks;
  }

  template <class ScalarT, class IdxT>
  int Exst1Pti<ScalarT, IdxT>::tagDifferentiable()
  {
    this->tag_[VC]   = Tr_ > 0;
    this->tag_[LL]   = Tb_ > 0;
    this->tag_[VA]   = Ta_ > 0;
    this->tag_[WASH] = Tf_ > 0;
    return 0;
  }

  template <class ScalarT, class IdxT>
  [[gnu::always_inline]] int Exst1Pti<ScalarT, IdxT>::evaluateInternalResidual(const ScalarT* y, const ScalarT* yp, const ScalarT* wb, const ScalarT* w, ScalarT* f)
  {
    using T           = ScalarT;
    const auto signal = [w](Input input) -> const T&
    { return w[static_cast<size_t>(input)]; };

    const T vt       = std::sqrt(wb[0] * wb[0] + wb[1] * wb[1]);
    const T feedback = washout_gain_ * (y[VA] - y[WASH]);
    const T error    = Math::clamp(signal(Input::vref) + signal(Input::vs) + signal(Input::vuel) + signal(Input::voel) - y[VC] - feedback, Vimin_, Vimax_);
    const T lead     = y[LL] + lead_gain_ * (error - y[LL]);
    const T drop     = Kc_ * signal(Input::ifd);

    f[VC]   = vt - y[VC] - Tr_ * yp[VC];
    f[LL]   = error - y[LL] - Tb_ * yp[LL];
    f[VA]   = Ka_ * lead - y[VA] - Ta_ * yp[VA];
    f[WASH] = y[VA] - y[WASH] - Tf_ * yp[WASH];
    f[EFD]  = Math::clamp(y[VA], vt * Vrmin_ - drop, vt * Vrmax_ - drop) - y[EFD];
    return 0;
  }

  template <class ScalarT, class IdxT>
  int Exst1Pti<ScalarT, IdxT>::initialize()
  {
    if (!this->allocated_ || !this->verifyAndLog())
      return 1;
    auto* y = this->y_.getData();
    this->yp_.setToZero();
    const auto input = [&](Input key, RealT fallback = 0)
    {
      const auto& port = ports_.in[key];
      return port ? static_cast<RealT>(port.readSignal()) : fallback;
    };

    const RealT vt    = std::hypot(static_cast<RealT>(bus_->Vr()), static_cast<RealT>(bus_->Vi()));
    const RealT efd   = static_cast<RealT>(y[EFD]);
    const RealT drop  = Kc_ * input(Input::ifd);
    const auto  field = [&](RealT error)
    {
      return Math::clamp(Ka_ * Math::clamp(error, Vimin_, Vimax_), vt * Vrmin_ - drop, vt * Vrmax_ - drop);
    };
    const RealT error = invertControl(field, efd, Vimin_ - 1, Vimax_ + 1);
    y[VC]             = vt;
    y[LL]             = Math::clamp(error, Vimin_, Vimax_);
    y[VA] = y[WASH] = Ka_ * y[LL];
    vref_set_       = vt + error - input(Input::vs) - input(Input::vuel) - input(Input::voel);
    if (auto& port = ports_.in[Input::vref]; port)
      port.writeValue(ScalarT{vref_set_});
    this->y_.setDataUpdated();
    evaluateResidual();
    for (size_t i = 0; i < COUNT; ++i)
    {
      const RealT r = static_cast<RealT>(this->f_.getData()[i]);
      if (!std::isfinite(r) || std::abs(r) > 1e-7 * (1 + std::abs(static_cast<RealT>(y[i]))))
      {
        Log::error() << "Exst1Pti: infeasible initial operating point, row " << i << ", residual " << r << '\n';
        return 1;
      }
    }
    if constexpr (std::is_same_v<ScalarT, DependencyTracking::Variable>)
      this->initializeDependencyTrackingVariableNumbers();
    return 0;
  }

  template <class ScalarT, class IdxT>
  int Exst1Pti<ScalarT, IdxT>::setSignalIndices()
  {
    for (auto input : Utilities::enum_values<Input>())
    {
      const auto& port                              = ports_.in[input];
      this->ws_indices_[static_cast<size_t>(input)] = port ? port.signalVariableIndex() : INVALID_INDEX<IdxT>;
    }
    return 0;
  }

  template <class ScalarT, class IdxT>
  int Exst1Pti<ScalarT, IdxT>::evaluateResidual()
  {
    auto*      w      = this->ws_.getData();
    const auto gather = [&](Input input, RealT fallback = 0)
    {
      const auto& port              = ports_.in[input];
      w[static_cast<size_t>(input)] = port ? port.readSignal() : ScalarT{fallback};
    };
    gather(Input::ifd);
    gather(Input::vref, vref_set_);
    gather(Input::vs);
    gather(Input::vuel);
    gather(Input::voel);
    this->wb_.getData()[0] = bus_->Vr();
    this->wb_.getData()[1] = bus_->Vi();
    evaluateInternalResidual(this->y_.getData(), this->yp_.getData(), this->wb_.getData(), w, this->f_.getData());
    this->f_.setDataUpdated();
    return 0;
  }

} // namespace GridKit::PhasorDynamics::Exciter
