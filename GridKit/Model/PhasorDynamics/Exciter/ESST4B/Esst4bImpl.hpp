#pragma once
#include <numbers>
#include <numeric>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/ParameterReader.hpp>
#include <GridKit/Model/PhasorDynamics/ControlLimits.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

#include "Esst4b.hpp"

namespace GridKit::PhasorDynamics::Exciter
{

  template <class ScalarT, class IdxT>
  Esst4b<ScalarT, IdxT>::Esst4b(BusT* bus, const DataT& data)
    : monitor_(std::make_unique<MonitorT>(data)), bus_(bus)
  {
    this->size_ = COUNT;
    monitor_->set(DataT::MonitorableVariables::efd, [this]
                  { return this->y_.getData()[EFD]; });
    using Parameters = typename DataT::Parameters;
    const GridKit::Model::ParameterReader reader(data, "Esst4b");
    mva_              = reader.readReal(Parameters::mva).value_or(mva_);
    Tr_               = reader.readReal(Parameters::Tr).value_or(Tr_);
    Kpr_              = reader.readReal(Parameters::Kpr).value_or(Kpr_);
    Kir_              = reader.readReal(Parameters::Kir).value_or(Kir_);
    Vrmax_            = reader.readReal(Parameters::Vrmax).value_or(Vrmax_);
    Vrmin_            = reader.readReal(Parameters::Vrmin).value_or(Vrmin_);
    Ta_               = reader.readReal(Parameters::Ta).value_or(Ta_);
    Kpm_              = reader.readReal(Parameters::Kpm).value_or(Kpm_);
    Kim_              = reader.readReal(Parameters::Kim).value_or(Kim_);
    Vmmax_            = reader.readReal(Parameters::Vmmax).value_or(Vmmax_);
    Vmmin_            = reader.readReal(Parameters::Vmmin).value_or(Vmmin_);
    Kg_               = reader.readReal(Parameters::Kg).value_or(Kg_);
    Kp_               = reader.readReal(Parameters::Kp).value_or(Kp_);
    Ki_               = reader.readReal(Parameters::Ki).value_or(Ki_);
    Vbmax_            = reader.readReal(Parameters::Vbmax).value_or(Vbmax_);
    Kc_               = reader.readReal(Parameters::Kc).value_or(Kc_);
    Xl_               = reader.readReal(Parameters::Xl).value_or(Xl_);
    ThetaP_           = reader.readReal(Parameters::ThetaP).value_or(ThetaP_);
    Vgmax_            = reader.readReal(Parameters::Vgmax).value_or(Vgmax_);
    const RealT angle = ThetaP_ * std::numbers::pi / 180;
    kp_real_          = Kp_ * std::cos(angle);
    kp_imag_          = Kp_ * std::sin(angle);

    this->setComponentBase(MEGA<RealT> * mva_);
    voel_set_ = Vmmax_ + 1;
  }

  template <class ScalarT, class IdxT>
  Esst4b<ScalarT, IdxT>::~Esst4b() = default;

  template <class ScalarT, class IdxT>
  int Esst4b<ScalarT, IdxT>::allocate()
  {
    compensated_ = ports_.in[Input::vcomp].connected();
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
  GridKit::Model::ConfigurationChecks Esst4b<ScalarT, IdxT>::verify() const
  {
    GridKit::Model::ConfigurationChecks checks("Esst4b");
    checks.check(bus_ != nullptr, "terminal bus required");
    for (auto input : Utilities::enum_values<Input>())
    {
      if (ports_.in[input] && !ports_.in[input].linked())
        checks.fail(std::format("unlinked input {}", Utilities::enum_name(input)));
    }

    checks.check(mva_ > 0, "mva must be positive");
    checks.check(Tr_ >= 0 && Ta_ >= 0, "Tr and Ta must be nonnegative");
    checks.check(Vrmin_ <= Vrmax_ && Vmmin_ <= Vmmax_, "regulator limits must be ordered");
    checks.check(Kc_ >= 0 && Vbmax_ > 0, "Kc must be nonnegative and Vbmax positive");
    checks.check(Kc_ == 0 || this->ports_.in[Input::ifd], "field current required with commutation loading");
    checks.check((Ki_ == 0 && Xl_ == 0) || (this->ports_.in[Input::ir] && this->ports_.in[Input::ii]), "compound excitation requires both terminal current components");
    checks.check(this->ports_.out[Output::efd].connected(), "efd output required");
    return checks;
  }

  template <class ScalarT, class IdxT>
  int Esst4b<ScalarT, IdxT>::tagDifferentiable()
  {
    this->tag_[VC]  = Tr_ > 0;
    this->tag_[XVR] = true;
    this->tag_[VA]  = Ta_ > 0;
    this->tag_[XVM] = true;
    return 0;
  }

  template <class ScalarT, class IdxT>
  template <class T>
  [[gnu::always_inline]] T Esst4b<ScalarT, IdxT>::bridgeVoltage(T vr, T vi, T ir, T ii, T ifd) const
  {
    const RealT a = kp_real_, b = kp_imag_;
    const T     er     = a * vr - b * vi - Xl_ * b * ir - (Ki_ + Xl_ * a) * ii;
    const T     ei     = b * vr + a * vi + (Ki_ + Xl_ * a) * ir - Xl_ * b * ii;
    const T     square = er * er + ei * ei;
    // V_E F_EX, evaluated on its three loading ranges without dividing by V_E.
    const T     ve     = square / std::sqrt(std::max(square, T{1e-24}));
    const T     load   = Kc_ * ifd;
    const T     knee   = (std::numbers::sqrt3 / 4) * ve;
    const T     first  = positiveDifference(load, T{0}) - positiveDifference(load, knee);
    const T     middle = knee + positiveDifference(load, knee) - positiveDifference(load, 0.75 * ve);
    const T     last   = 0.75 * ve + positiveDifference(load, 0.75 * ve) - positiveDifference(load, ve);
    const T     arc    = 0.75 * square - middle * middle;
    const T     bridge = 0.25 * ve - first / std::numbers::sqrt3 + arc / std::sqrt(std::max(arc, T{1e-24}))
                     - std::numbers::sqrt3 * (last - 0.75 * ve);

    return Math::min(bridge, Vbmax_);
  }

  template <class ScalarT, class IdxT>
  [[gnu::always_inline]] int Esst4b<ScalarT, IdxT>::evaluateInternalResidual(const ScalarT* y, const ScalarT* yp, const ScalarT* wb, const ScalarT* w, ScalarT* f)
  {
    using T           = ScalarT;
    const auto signal = [w](Input input) -> const T&
    { return w[static_cast<size_t>(input)]; };

    const T bridge       = bridgeVoltage(wb[0], wb[1], this->toComponentBase(signal(Input::ir)), this->toComponentBase(signal(Input::ii)), signal(Input::ifd));
    const T error        = signal(Input::vref) + signal(Input::vs) + signal(Input::vuel) - y[VC];
    const T vr_unlimited = Kpr_ * error + y[XVR];
    const T field_error  = y[VA] - Math::min(Kg_ * y[EFD], Vgmax_);
    const T vm           = Kpm_ * field_error + y[XVM];

    f[VC]  = (compensated_ * signal(Input::vcomp) + (1 - compensated_) * std::sqrt(wb[0] * wb[0] + wb[1] * wb[1])) - y[VC] - Tr_ * yp[VC];
    f[XVR] = Math::antiwindup(vr_unlimited, Kir_ * error, Vrmin_, Vrmax_) - yp[XVR];
    f[VA]  = Math::clamp(vr_unlimited, Vrmin_, Vrmax_) - y[VA] - Ta_ * yp[VA];
    f[XVM] = Math::antiwindup(vm, Kim_ * field_error, Vmmin_, Vmmax_) - yp[XVM];
    f[EFD] = bridge * Math::min(Math::clamp(vm, Vmmin_, Vmmax_), signal(Input::voel)) - y[EFD];
    return 0;
  }

  template <class ScalarT, class IdxT>
  int Esst4b<ScalarT, IdxT>::initialize()
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

    const RealT vr = static_cast<RealT>(bus_->Vr()), vi = static_cast<RealT>(bus_->Vi());
    const RealT vc     = compensated_ != 0 ? input(Input::vcomp) : std::hypot(vr, vi);
    const RealT efd    = static_cast<RealT>(y[EFD]);
    const RealT bridge = bridgeVoltage(vr, vi, this->toComponentBase(input(Input::ir)), this->toComponentBase(input(Input::ii)), input(Input::ifd));
    if (bridge <= 0 && efd != 0)
      return 1;
    const RealT va = Math::min(Kg_ * efd, Vgmax_);
    y[VC]          = vc;
    y[VA]          = va;
    y[XVR]         = invertControl([&](RealT x)
                           { return Math::clamp(x, Vrmin_, Vrmax_); },
                           va,
                           Vrmin_ - 1,
                           Vrmax_ + 1);
    y[XVM]         = bridge <= 0 ? 0 : invertControl([&](RealT x)
                                             { return Math::min(Math::clamp(x, Vmmin_, Vmmax_), input(Input::voel, voel_set_)); },
                                             efd / bridge,
                                             Vmmin_ - 1,
                                             Vmmax_ + 1);
    vref_set_      = vc - input(Input::vs) - input(Input::vuel);
    if (auto& port = ports_.in[Input::vref]; port)
      port.writeValue(ScalarT{vref_set_});
    this->y_.setDataUpdated();
    evaluateResidual();
    for (size_t i = 0; i < COUNT; ++i)
    {
      const RealT r = static_cast<RealT>(this->f_.getData()[i]);
      if (!std::isfinite(r) || std::abs(r) > 1e-7 * (1 + std::abs(static_cast<RealT>(y[i]))))
      {
        Log::error() << "Esst4b: infeasible initial operating point, row " << i << ", residual " << r << '\n';
        return 1;
      }
    }
    if constexpr (std::is_same_v<ScalarT, DependencyTracking::Variable>)
      this->initializeDependencyTrackingVariableNumbers();
    return 0;
  }

  template <class ScalarT, class IdxT>
  int Esst4b<ScalarT, IdxT>::setSignalIndices()
  {
    for (auto input : Utilities::enum_values<Input>())
    {
      const auto& port                              = ports_.in[input];
      this->ws_indices_[static_cast<size_t>(input)] = port ? port.signalVariableIndex() : INVALID_INDEX<IdxT>;
    }
    return 0;
  }

  template <class ScalarT, class IdxT>
  int Esst4b<ScalarT, IdxT>::evaluateResidual()
  {
    auto*      w      = this->ws_.getData();
    const auto gather = [&](Input input, RealT fallback = 0)
    {
      const auto& port              = ports_.in[input];
      w[static_cast<size_t>(input)] = port ? port.readSignal() : ScalarT{fallback};
    };
    gather(Input::ifd);
    gather(Input::ir);
    gather(Input::ii);
    gather(Input::vref, vref_set_);
    gather(Input::vs);
    gather(Input::vuel);
    gather(Input::voel, voel_set_);
    gather(Input::vcomp);
    this->wb_.getData()[0] = bus_->Vr();
    this->wb_.getData()[1] = bus_->Vi();
    evaluateInternalResidual(this->y_.getData(), this->yp_.getData(), this->wb_.getData(), w, this->f_.getData());
    this->f_.setDataUpdated();
    return 0;
  }

} // namespace GridKit::PhasorDynamics::Exciter
