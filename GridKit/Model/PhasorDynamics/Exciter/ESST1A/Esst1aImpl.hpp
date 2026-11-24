#pragma once
#include <numeric>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/ParameterReader.hpp>
#include <GridKit/Model/PhasorDynamics/ControlLimits.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

#include "Esst1a.hpp"

namespace GridKit::PhasorDynamics::Exciter
{

  template <class ScalarT, class IdxT>
  Esst1a<ScalarT, IdxT>::Esst1a(BusT* bus, const DataT& data)
    : monitor_(std::make_unique<MonitorT>(data)), bus_(bus)
  {
    this->size_ = COUNT;
    monitor_->set(DataT::MonitorableVariables::efd, [this]
                  { return this->y_.getData()[EFD]; });
    using Parameters = typename DataT::Parameters;
    const GridKit::Model::ParameterReader reader(data, "Esst1a");
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
    UEL_   = reader.readInteger(Parameters::UEL).value_or(UEL_);
    VOS_   = reader.readInteger(Parameters::VOS).value_or(VOS_);
    Tc1_   = reader.readReal(Parameters::Tc1).value_or(Tc1_);
    Tb1_   = reader.readReal(Parameters::Tb1).value_or(Tb1_);
    Vamax_ = reader.readReal(Parameters::Vamax).value_or(Vamax_);
    Vamin_ = reader.readReal(Parameters::Vamin).value_or(Vamin_);
    Klr_   = reader.readReal(Parameters::Klr).value_or(Klr_);
    Ilr_   = reader.readReal(Parameters::Ilr).value_or(Ilr_);

    washout_gain_     = Tf_ > 0 ? Kf_ / Tf_ : 0;
    lead_gain_        = Tb_ > 0 ? Tc_ / Tb_ : 0;
    second_lead_gain_ = Tb1_ > 0 ? Tc1_ / Tb1_ : 0;

    voel_set_ = std::max({1.0, Vamax_, Vrmax_}) + 1;
  }

  template <class ScalarT, class IdxT>
  Esst1a<ScalarT, IdxT>::~Esst1a() = default;

  template <class ScalarT, class IdxT>
  int Esst1a<ScalarT, IdxT>::allocate()
  {
    gate_input_  = ports_.in[Input::vuel].connected() && UEL_ == 2;
    gate_output_ = ports_.in[Input::vuel].connected() && UEL_ == 3;
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
  GridKit::Model::ConfigurationChecks Esst1a<ScalarT, IdxT>::verify() const
  {
    GridKit::Model::ConfigurationChecks checks("Esst1a");
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
    checks.check(Tc1_ >= 0 && Tb1_ >= 0 && (Tb1_ > 0 || Tc1_ == 0), "Tc1 requires a positive Tb1");
    checks.check(Vamin_ <= Vamax_, "regulator limits must be ordered");
    checks.check(static_cast<size_t>(UEL_) <= 3 && (VOS_ == 1 || VOS_ == 2), "invalid UEL or VOS input location");
    checks.check((Kc_ == 0 && Klr_ == 0) || this->ports_.in[Input::ifd], "field current required for loading or limiting");
    checks.check(this->ports_.out[Output::efd].connected(), "efd output required");
    return checks;
  }

  template <class ScalarT, class IdxT>
  int Esst1a<ScalarT, IdxT>::tagDifferentiable()
  {
    this->tag_[VC]   = Tr_ > 0;
    this->tag_[LL]   = Tb_ > 0;
    this->tag_[LL1]  = Tb1_ > 0;
    this->tag_[VA]   = Ta_ > 0;
    this->tag_[WASH] = Tf_ > 0;
    return 0;
  }

  template <class ScalarT, class IdxT>
  [[gnu::always_inline]] int Esst1a<ScalarT, IdxT>::evaluateInternalResidual(const ScalarT* y, const ScalarT* yp, const ScalarT* wb, const ScalarT* w, ScalarT* f)
  {
    using T           = ScalarT;
    const auto signal = [w](Input input) -> const T&
    { return w[static_cast<size_t>(input)]; };

    const T vt       = std::sqrt(wb[0] * wb[0] + wb[1] * wb[1]);
    const T ifd      = signal(Input::ifd);
    const T limited  = y[VA] + (VOS_ == 2) * signal(Input::vs) - Klr_ * Math::ramp(ifd - Ilr_);
    const T gated    = limited + gate_output_ * Math::ramp(signal(Input::vuel) - limited);
    f[REG]           = Math::min(gated, signal(Input::voel)) - y[REG];
    const T feedback = washout_gain_ * (y[REG] - y[WASH]);
    const T error    = signal(Input::vref) - y[VC] - feedback + (VOS_ == 1) * signal(Input::vs)
                    + (UEL_ == 1) * signal(Input::vuel);
    const T vi = Math::clamp(error, Vimin_, Vimax_);
    f[DRIVE]   = vi + gate_input_ * Math::ramp(signal(Input::vuel) - vi) - y[DRIVE];
    f[STAGE]   = y[LL] + lead_gain_ * (y[DRIVE] - y[LL]) - y[STAGE];
    f[DRIVER]  = Ka_ * (y[LL1] + second_lead_gain_ * (y[STAGE] - y[LL1])) - y[DRIVER];

    f[VC]  = vt - y[VC] - Tr_ * yp[VC];
    f[LL]  = y[DRIVE] - y[LL] - Tb_ * yp[LL];
    f[LL1] = y[STAGE] - y[LL1] - Tb1_ * yp[LL1];
    f[VA]  = (Ta_ > 0) * Math::antiwindup(y[VA], y[DRIVER] - y[VA], Vamin_, Vamax_)
            + (Ta_ == 0) * (Math::clamp(y[DRIVER], Vamin_, Vamax_) - y[VA]) - Ta_ * yp[VA];
    f[WASH] = y[REG] - y[WASH] - Tf_ * yp[WASH];
    f[EFD]  = Math::max(vt * Vrmin_, Math::min(y[REG], vt * Vrmax_ - Kc_ * ifd)) - y[EFD];
    return 0;
  }

  template <class ScalarT, class IdxT>
  int Esst1a<ScalarT, IdxT>::initialize()
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

    const RealT vt  = std::hypot(static_cast<RealT>(bus_->Vr()), static_cast<RealT>(bus_->Vi()));
    const RealT ifd = input(Input::ifd), vs = input(Input::vs), vuel = input(Input::vuel), voel = input(Input::voel, voel_set_);
    const RealT efd      = static_cast<RealT>(y[EFD]);
    const auto  regulate = [&](RealT error)
    {
      const RealT vi      = Math::clamp(error, Vimin_, Vimax_);
      const RealT drive   = vi + gate_input_ * Math::ramp(vuel - vi);
      const RealT driver  = Ka_ * drive;
      const RealT va      = Ta_ > 0 ? std::clamp(driver, Vamin_, Vamax_) : Math::clamp(driver, Vamin_, Vamax_);
      const RealT limited = va + (VOS_ == 2) * vs - Klr_ * Math::ramp(ifd - Ilr_);
      const RealT reg     = Math::min(limited + gate_output_ * Math::ramp(vuel - limited), voel);
      return std::array{drive, driver, va, reg};
    };
    const auto field = [&](RealT error)
    {
      return Math::max(vt * Vrmin_, Math::min(regulate(error)[3], vt * Vrmax_ - Kc_ * ifd));
    };
    const RealT error                   = invertControl(field, efd, Vimin_ - 1, Vimax_ + 1);
    const auto [drive, driver, va, reg] = regulate(error);
    y[VC]                               = vt;
    y[LL] = y[LL1] = y[STAGE] = y[DRIVE] = drive;
    y[DRIVER]                            = driver;
    y[VA]                                = va;
    y[WASH] = y[REG] = reg;
    vref_set_        = vt + error - (VOS_ == 1) * vs - (UEL_ == 1) * vuel;
    if (auto& port = ports_.in[Input::vref]; port)
      port.writeValue(ScalarT{vref_set_});
    this->y_.setDataUpdated();
    evaluateResidual();
    for (size_t i = 0; i < COUNT; ++i)
    {
      const RealT r = static_cast<RealT>(this->f_.getData()[i]);
      if (!std::isfinite(r) || std::abs(r) > 1e-7 * (1 + std::abs(static_cast<RealT>(y[i]))))
      {
        Log::error() << "Esst1a: infeasible initial operating point, row " << i << ", residual " << r << '\n';
        return 1;
      }
    }
    if constexpr (std::is_same_v<ScalarT, DependencyTracking::Variable>)
      this->initializeDependencyTrackingVariableNumbers();
    return 0;
  }

  template <class ScalarT, class IdxT>
  int Esst1a<ScalarT, IdxT>::setSignalIndices()
  {
    for (auto input : Utilities::enum_values<Input>())
    {
      const auto& port                              = ports_.in[input];
      this->ws_indices_[static_cast<size_t>(input)] = port ? port.signalVariableIndex() : INVALID_INDEX<IdxT>;
    }
    return 0;
  }

  template <class ScalarT, class IdxT>
  int Esst1a<ScalarT, IdxT>::evaluateResidual()
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
    gather(Input::voel, voel_set_);
    this->wb_.getData()[0] = bus_->Vr();
    this->wb_.getData()[1] = bus_->Vi();
    evaluateInternalResidual(this->y_.getData(), this->yp_.getData(), this->wb_.getData(), w, this->f_.getData());
    this->f_.setDataUpdated();
    return 0;
  }

} // namespace GridKit::PhasorDynamics::Exciter
