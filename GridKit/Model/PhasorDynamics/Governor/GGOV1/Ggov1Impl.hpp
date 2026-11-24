#pragma once
#include <numeric>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/ParameterReader.hpp>
#include <GridKit/Model/PhasorDynamics/ControlLimits.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

#include "Ggov1.hpp"

namespace GridKit::PhasorDynamics::Governor
{

  template <class ScalarT, class IdxT>
  Ggov1<ScalarT, IdxT>::Ggov1(const DataT& data)
    : monitor_(std::make_unique<MonitorT>(data))
  {
    this->size_ = COUNT;
    monitor_->set(DataT::MonitorableVariables::pmech, [this]
                  { return this->y_.getData()[PM]; });
    monitor_->set(DataT::MonitorableVariables::valve, [this]
                  { return this->y_.getData()[VALVE]; });
    using Parameters = typename DataT::Parameters;
    const GridKit::Model::ParameterReader reader(data, "Ggov1");
    mva_   = reader.readReal(Parameters::mva).value_or(mva_);
    Trate_ = reader.readReal(Parameters::Trate).value_or(Trate_);
    if (auto value = reader.readInteger(Parameters::Rselect))
    {
      // JSON integers use IdxT, including the two negative feedback selectors.
      Rselect_ = static_cast<std::make_signed_t<IdxT>>(*value);
    }
    Flag_              = reader.readBool(Parameters::Flag).value_or(Flag_);
    R_                 = reader.readReal(Parameters::R).value_or(R_);
    Tpelec_            = reader.readReal(Parameters::Tpelec).value_or(Tpelec_);
    Maxerr_            = reader.readReal(Parameters::Maxerr).value_or(Maxerr_);
    Minerr_            = reader.readReal(Parameters::Minerr).value_or(Minerr_);
    Kpgov_             = reader.readReal(Parameters::Kpgov).value_or(Kpgov_);
    Kigov_             = reader.readReal(Parameters::Kigov).value_or(Kigov_);
    Kdgov_             = reader.readReal(Parameters::Kdgov).value_or(Kdgov_);
    Tdgov_             = reader.readReal(Parameters::Tdgov).value_or(Tdgov_);
    Vmax_              = reader.readReal(Parameters::Vmax).value_or(Vmax_);
    Vmin_              = reader.readReal(Parameters::Vmin).value_or(Vmin_);
    Tact_              = reader.readReal(Parameters::Tact).value_or(Tact_);
    Kturb_             = reader.readReal(Parameters::Kturb).value_or(Kturb_);
    Wfnl_              = reader.readReal(Parameters::Wfnl).value_or(Wfnl_);
    Tb_                = reader.readReal(Parameters::Tb).value_or(Tb_);
    Tc_                = reader.readReal(Parameters::Tc).value_or(Tc_);
    Teng_              = reader.readReal(Parameters::Teng).value_or(Teng_);
    Tfload_            = reader.readReal(Parameters::Tfload).value_or(Tfload_);
    Kpload_            = reader.readReal(Parameters::Kpload).value_or(Kpload_);
    Kiload_            = reader.readReal(Parameters::Kiload).value_or(Kiload_);
    Ldref_             = reader.readReal(Parameters::Ldref).value_or(Ldref_);
    Dm_                = reader.readReal(Parameters::Dm).value_or(Dm_);
    Ropen_             = reader.readReal(Parameters::Ropen).value_or(Ropen_);
    Rclose_            = reader.readReal(Parameters::Rclose).value_or(Rclose_);
    Kimw_              = reader.readReal(Parameters::Kimw).value_or(Kimw_);
    Aset_              = reader.readReal(Parameters::Aset).value_or(Aset_);
    Ka_                = reader.readReal(Parameters::Ka).value_or(Ka_);
    Ta_                = reader.readReal(Parameters::Ta).value_or(Ta_);
    Db_                = reader.readReal(Parameters::Db).value_or(Db_);
    Tsa_               = reader.readReal(Parameters::Tsa).value_or(Tsa_);
    Tsb_               = reader.readReal(Parameters::Tsb).value_or(Tsb_);
    Tfsr_              = reader.readReal(Parameters::Tfsr).value_or(Tfsr_);
    Fsrtmax_           = reader.readReal(Parameters::Fsrtmax).value_or(Fsrtmax_);
    governor_tracking_ = Kpgov_ > 0 ? Kigov_ / Kpgov_ : 0;
    load_tracking_     = Kpload_ > 0 ? Kiload_ / Kpload_ : 0;
    temperature_gain_  = Tsb_ > 0 ? Tsa_ / Tsb_ : 0;
    turbine_gain_      = Tb_ > 0 ? Tc_ / Tb_ : 0;
    acceleration_gain_ = Ta_ > 0 ? Ka_ / Ta_ : 0;
    derivative_gain_   = Tdgov_ > 0 ? Kdgov_ / Tdgov_ : 0;

    this->setComponentBase(MEGA<RealT> * (Trate_ > 0 ? Trate_ : mva_));
  }

  template <class ScalarT, class IdxT>
  Ggov1<ScalarT, IdxT>::~Ggov1() = default;

  template <class ScalarT, class IdxT>
  int Ggov1<ScalarT, IdxT>::allocate()
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
    if (auto& port = ports_.out[Output::pmech]; port)
      port.link(this->y_.getData() + PM, &this->getVariableIndex(PM));
    this->allocated_ = true;
    return 0;
  }

  template <class ScalarT, class IdxT>
  GridKit::Model::ConfigurationChecks Ggov1<ScalarT, IdxT>::verify() const
  {
    GridKit::Model::ConfigurationChecks checks("Ggov1");
    for (auto input : Utilities::enum_values<Input>())
    {
      if (ports_.in[input] && !ports_.in[input].linked())
        checks.fail(std::format("unlinked input {}", Utilities::enum_name(input)));
    }

    checks.check((Trate_ > 0 || mva_ > 0) && Trate_ >= 0, "Trate must be nonnegative and requires positive mva when zero");
    checks.check(-2 <= Rselect_ && Rselect_ <= 1 && R_ >= 0, "invalid droop selection or gain");
    checks.check(Kpgov_ >= 0 && Kigov_ >= 0 && Kpload_ >= 0 && Kiload_ >= 0, "PI gains must be nonnegative");
    checks.check(Tpelec_ >= 0 && Tdgov_ >= 0 && (Tdgov_ > 0 || Kdgov_ == 0), "invalid measurement or derivative time");
    checks.check(Tact_ > 0 && Kturb_ > 0 && Tfload_ >= 0, "invalid actuator, turbine or load limiter");
    checks.check(Tb_ >= 0 && Tc_ >= 0 && (Tb_ > 0 || Tc_ == 0), "Tc requires a positive Tb");
    checks.check(Tsb_ >= 0 && Tsa_ >= 0 && (Tsb_ > 0 || Tsa_ == 0), "Tsa requires a positive Tsb");
    checks.check(Teng_ == 0 && Tfsr_ > 0, "GGOV1 requires zero Teng and positive fuel memory Tfsr");
    checks.check(Ta_ >= 0 && (Ka_ == 0 || Ta_ > 0), "acceleration control requires a positive Ta");
    checks.check(Minerr_ <= Maxerr_ && Vmin_ <= Vmax_ && Rclose_ < 0 && Ropen_ > 0, "invalid governor limits");
    checks.check(this->ports_.in[Input::pelec].connected(), "pelec input required");
    checks.check(this->ports_.out[Output::pmech].connected(), "pmech output required");
    return checks;
  }

  template <class ScalarT, class IdxT>
  int Ggov1<ScalarT, IdxT>::tagDifferentiable()
  {
    this->tag_[PE]     = Tpelec_ > 0;
    this->tag_[XD]     = Tdgov_ > 0;
    this->tag_[XI]     = true;
    this->tag_[VALVE]  = true;
    this->tag_[TURB]   = Tb_ > 0;
    this->tag_[TEMP]   = Tsb_ > 0;
    this->tag_[LOAD]   = Tfload_ > 0;
    this->tag_[XL]     = true;
    this->tag_[XMW]    = true;
    this->tag_[SPEED]  = Ta_ > 0;
    this->tag_[MEMORY] = true;
    return 0;
  }

  template <class ScalarT, class IdxT>
  [[gnu::always_inline]] int Ggov1<ScalarT, IdxT>::evaluateInternalResidual(const ScalarT* y, const ScalarT* yp, const ScalarT*, const ScalarT* w, ScalarT* f)
  {
    using T           = ScalarT;
    const auto signal = [w](Input input) -> const T&
    { return w[static_cast<size_t>(input)]; };

    const T speed      = signal(Input::speed);
    const T feedback   = (Rselect_ == 1) * y[PE] + (Rselect_ == -1) * y[VALVE] + (Rselect_ == -2) * y[FSR];
    f[ERROR]           = Math::clamp(Math::deadband2(this->toComponentBase(signal(Input::pref)) + y[XMW] - speed - R_ * feedback, -Db_, Db_), Minerr_, Maxerr_) - y[ERROR];
    const T derivative = derivative_gain_ * (y[ERROR] - y[XD]);
    f[FSRN]            = Kpgov_ * y[ERROR] + y[XI] + derivative - y[FSRN];
    const T fuel       = y[VALVE] * (ONE<RealT> + Flag_ * speed);
    const T exhaust    = fuel * std::exp(std::min(Dm_, RealT{0}) * std::log(Math::max((ONE<RealT> + speed), RealT{0.01})));
    const T sensed     = y[TEMP] + temperature_gain_ * (exhaust - y[TEMP]);
    const T load_error = (Ldref_ + this->toComponentBase(signal(Input::paux))) / Kturb_ + Wfnl_ - y[LOAD];
    f[FSRT]            = Kpload_ * load_error + y[XL] - y[FSRT];
    const T turbine    = y[TURB] + turbine_gain_ * (fuel - y[TURB]);

    f[PE]               = this->toComponentBase(signal(Input::pelec)) - y[PE] - Tpelec_ * yp[PE];
    f[XD]               = y[ERROR] - y[XD] - Tdgov_ * yp[XD];
    f[XI]               = governor_tracking_ * (y[FSR] - y[XI]) + (Kpgov_ == 0) * Kigov_ * y[ERROR] - yp[XI];
    f[VALVE]            = Math::clamp((y[FSR] - y[VALVE]) / Tact_, Rclose_, Ropen_) - yp[VALVE];
    f[TURB]             = fuel - y[TURB] - Tb_ * yp[TURB];
    f[TEMP]             = exhaust - y[TEMP] - Tsb_ * yp[TEMP];
    f[LOAD]             = sensed - y[LOAD] - Tfload_ * yp[LOAD];
    f[XL]               = load_tracking_ * (y[FSR] - y[XL]) + (Kpload_ == 0) * Kiload_ * load_error - yp[XL];
    f[XMW]              = Math::antiwindup(y[XMW], Kimw_ * (this->toComponentBase(signal(Input::pmwset)) - y[PE]), -1.1 * R_, 1.1 * R_) - yp[XMW];
    f[SPEED]            = speed - y[SPEED] - Ta_ * yp[SPEED];
    // Exact fuel selection preserves the operating point at coincident orders and valve limits.
    const auto positive = [](T x)
    { return HALF<RealT> * (x + std::abs(x)); };
    const auto minimum = [&](T a, T b)
    { return a - positive(a - b); };
    const auto limited = [&](T x)
    { return Vmin_ + positive(x - Vmin_) - positive(x - Vmax_); };
    f[DEMAND]     = limited(minimum(y[FSRN], minimum(y[FSRT], T{Fsrtmax_}))) - y[DEMAND];
    // The acceleration controller compares demand with the previous fuel command.
    const T accel = Ka_ * Aset_ - acceleration_gain_ * (speed - y[SPEED]);
    const T order = y[DEMAND] - (Ka_ > 0) * positive(y[DEMAND] - y[MEMORY] - Tfsr_ * accel);
    f[FSR]        = limited(order) - y[FSR];
    f[MEMORY]     = y[FSR] - y[MEMORY] - Tfsr_ * yp[MEMORY];
    f[PM]         = Kturb_ * (turbine - Wfnl_) - std::max(Dm_, RealT{0}) * speed - this->toComponentBase(y[PM]);
    return 0;
  }

  template <class ScalarT, class IdxT>
  int Ggov1<ScalarT, IdxT>::initialize()
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

    const RealT speed         = input(Input::speed);
    const RealT pe            = this->toComponentBase(input(Input::pelec));
    const RealT mechanical    = this->toComponentBase(static_cast<RealT>(y[PM]));
    const RealT fuel          = (mechanical + std::max(Dm_, RealT{0}) * speed) / Kturb_ + Wfnl_;
    const RealT valve         = fuel / (1 + Flag_ * speed);
    const RealT exhaust       = fuel * std::pow(std::max(1 + speed, RealT{0.01}), std::min(Dm_, RealT{0}));
    const auto  actuator_rate = [&](RealT x)
    { return Math::clamp(x, Rclose_, Ropen_); };
    const RealT fsr            = valve + Tact_ * invertControl(actuator_rate, RealT{0}, Rclose_ - 1, Ropen_ + 1);
    const RealT feedback       = (Rselect_ == 1) * pe + (Rselect_ == -1) * valve + (Rselect_ == -2) * fsr;
    const auto  governor_error = [&](RealT x)
    { return Math::clamp(Math::deadband2(x, -Db_, Db_), Minerr_, Maxerr_); };
    const RealT error_offset = invertControl(governor_error, RealT{0}, Minerr_ - Db_ - 1, Maxerr_ + Db_ + 1);
    y[PE]                    = pe;
    y[VALVE]                 = valve;
    y[FSR] = y[FSRN] = y[DEMAND] = y[MEMORY] = fsr;
    y[XI]                                    = fsr;
    y[TEMP] = y[LOAD] = exhaust;
    y[TURB]           = fuel;
    y[SPEED]          = speed;
    y[ERROR] = y[XD] = y[XMW] = 0;
    const RealT load_error    = (Ldref_ + this->toComponentBase(input(Input::paux))) / Kturb_ + Wfnl_ - exhaust;
    y[XL]                     = fsr;
    y[FSRT]                   = fsr + Kpload_ * load_error;
    pref_set_                 = this->toSystemBase(speed + R_ * feedback + error_offset);
    if (auto& port = ports_.in[Input::pref]; port)
      port.writeValue(ScalarT{pref_set_});
    pmwset_set_ = this->toSystemBase(pe);
    if (auto& port = ports_.in[Input::pmwset]; port)
      port.writeValue(ScalarT{pmwset_set_});
    this->y_.setDataUpdated();
    evaluateResidual();
    for (size_t i = 0; i < COUNT; ++i)
    {
      const RealT r = static_cast<RealT>(this->f_.getData()[i]);
      if (!std::isfinite(r) || std::abs(r) > 1e-7 * (1 + std::abs(static_cast<RealT>(y[i]))))
      {
        Log::error() << "Ggov1: infeasible initial operating point, row " << i << ", residual " << r << '\n';
        return 1;
      }
    }
    if constexpr (std::is_same_v<ScalarT, DependencyTracking::Variable>)
      this->initializeDependencyTrackingVariableNumbers();
    return 0;
  }

  template <class ScalarT, class IdxT>
  int Ggov1<ScalarT, IdxT>::setSignalIndices()
  {
    for (auto input : Utilities::enum_values<Input>())
    {
      const auto& port                              = ports_.in[input];
      this->ws_indices_[static_cast<size_t>(input)] = port ? port.signalVariableIndex() : INVALID_INDEX<IdxT>;
    }
    return 0;
  }

  template <class ScalarT, class IdxT>
  int Ggov1<ScalarT, IdxT>::evaluateResidual()
  {
    auto*      w      = this->ws_.getData();
    const auto gather = [&](Input input, RealT fallback = 0)
    {
      const auto& port              = ports_.in[input];
      w[static_cast<size_t>(input)] = port ? port.readSignal() : ScalarT{fallback};
    };
    gather(Input::speed);
    gather(Input::pelec);
    gather(Input::pref, pref_set_);
    gather(Input::paux);
    gather(Input::pmwset, pmwset_set_);
    evaluateInternalResidual(this->y_.getData(), this->yp_.getData(), nullptr, w, this->f_.getData());
    this->f_.setDataUpdated();
    return 0;
  }

} // namespace GridKit::PhasorDynamics::Governor
