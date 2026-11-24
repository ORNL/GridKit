#pragma once
#include <numeric>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/CurrentCircle.hpp>
#include <GridKit/Model/ParameterReader.hpp>
#include <GridKit/Model/PhasorDynamics/ControlLimits.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

#include "Reeca.hpp"

namespace GridKit::PhasorDynamics::Controller
{

  template <class ScalarT, class IdxT>
  Reeca<ScalarT, IdxT>::Reeca(BusT* bus, const DataT& data)
    : monitor_(std::make_unique<MonitorT>(data)), bus_(bus)
  {
    this->size_ = COUNT;
    monitor_->set(DataT::MonitorableVariables::iqcmd, [this]
                  { return this->y_.getData()[IQCMD]; });
    monitor_->set(DataT::MonitorableVariables::ipcmd, [this]
                  { return this->y_.getData()[IPCMD]; });
    monitor_->set(DataT::MonitorableVariables::vt, [this]
                  { return this->y_.getData()[VTERM]; });
    monitor_->set(DataT::MonitorableVariables::ipmax, [this]
                  { return this->y_.getData()[IPMAX]; });
    monitor_->set(DataT::MonitorableVariables::iqinj, [this]
                  { return this->y_.getData()[IQV]; });
    using Parameters = typename DataT::Parameters;
    const GridKit::Model::ParameterReader reader(data, "Reeca");
    mva_    = reader.readReal(Parameters::mva).value_or(mva_);
    PfFlag_ = reader.readBool(Parameters::PfFlag).value_or(PfFlag_);
    VFlag_  = reader.readBool(Parameters::VFlag).value_or(VFlag_);
    QFlag_  = reader.readBool(Parameters::QFlag).value_or(QFlag_);
    Pqflag_ = reader.readBool(Parameters::Pqflag).value_or(Pqflag_);
    Pflag_  = reader.readBool(Parameters::Pflag).value_or(Pflag_);
    Trv_    = reader.readReal(Parameters::Trv).value_or(Trv_);
    Tp_     = reader.readReal(Parameters::Tp).value_or(Tp_);
    Vref0_  = reader.readReal(Parameters::Vref0).value_or(Vref0_);
    Vref1_  = reader.readReal(Parameters::Vref1).value_or(Vref1_);
    Vdip_   = reader.readReal(Parameters::Vdip).value_or(Vdip_);
    Vup_    = reader.readReal(Parameters::Vup).value_or(Vup_);
    dbd1_   = reader.readReal(Parameters::dbd1).value_or(dbd1_);
    dbd2_   = reader.readReal(Parameters::dbd2).value_or(dbd2_);
    kqv_    = reader.readReal(Parameters::kqv).value_or(kqv_);
    Iql1_   = reader.readReal(Parameters::Iql1).value_or(Iql1_);
    Iqh1_   = reader.readReal(Parameters::Iqh1).value_or(Iqh1_);
    Iqfrz_  = reader.readReal(Parameters::Iqfrz).value_or(Iqfrz_);
    Thld_   = reader.readReal(Parameters::Thld).value_or(Thld_);
    Thld2_  = reader.readReal(Parameters::Thld2).value_or(Thld2_);
    Qmax_   = reader.readReal(Parameters::Qmax).value_or(Qmax_);
    Qmin_   = reader.readReal(Parameters::Qmin).value_or(Qmin_);
    Kqp_    = reader.readReal(Parameters::Kqp).value_or(Kqp_);
    Kqi_    = reader.readReal(Parameters::Kqi).value_or(Kqi_);
    Vmax_   = reader.readReal(Parameters::Vmax).value_or(Vmax_);
    Vmin_   = reader.readReal(Parameters::Vmin).value_or(Vmin_);
    Kvp_    = reader.readReal(Parameters::Kvp).value_or(Kvp_);
    Kvi_    = reader.readReal(Parameters::Kvi).value_or(Kvi_);
    Tiq_    = reader.readReal(Parameters::Tiq).value_or(Tiq_);
    Tpord_  = reader.readReal(Parameters::Tpord).value_or(Tpord_);
    dPmax_  = reader.readReal(Parameters::dPmax).value_or(dPmax_);
    dPmin_  = reader.readReal(Parameters::dPmin).value_or(dPmin_);
    Pmax_   = reader.readReal(Parameters::Pmax).value_or(Pmax_);
    Pmin_   = reader.readReal(Parameters::Pmin).value_or(Pmin_);
    Imax_   = reader.readReal(Parameters::Imax).value_or(Imax_);
    Vq_[0]  = reader.readReal(Parameters::Vq1).value_or(0);
    Vq_[1]  = reader.readReal(Parameters::Vq2).value_or(0);
    Vq_[2]  = reader.readReal(Parameters::Vq3).value_or(0);
    Vq_[3]  = reader.readReal(Parameters::Vq4).value_or(0);
    Iq_[0]  = reader.readReal(Parameters::Iq1).value_or(0);
    Iq_[1]  = reader.readReal(Parameters::Iq2).value_or(0);
    Iq_[2]  = reader.readReal(Parameters::Iq3).value_or(0);
    Iq_[3]  = reader.readReal(Parameters::Iq4).value_or(0);
    Vp_[0]  = reader.readReal(Parameters::Vp1).value_or(0);
    Vp_[1]  = reader.readReal(Parameters::Vp2).value_or(0);
    Vp_[2]  = reader.readReal(Parameters::Vp3).value_or(0);
    Vp_[3]  = reader.readReal(Parameters::Vp4).value_or(0);
    Ip_[0]  = reader.readReal(Parameters::Ip1).value_or(0);
    Ip_[1]  = reader.readReal(Parameters::Ip2).value_or(0);
    Ip_[2]  = reader.readReal(Parameters::Ip3).value_or(0);
    Ip_[3]  = reader.readReal(Parameters::Ip4).value_or(0);

    for (size_t k = 0; k < slope_q_.size(); ++k)
    {
      const auto slope = [k](const auto& v, const auto& i)
      {
        return v[k + 1] > v[k] ? (i[k + 1] - i[k]) / (v[k + 1] - v[k]) : RealT{0};
      };
      slope_q_[k] = slope(Vq_, Iq_);
      slope_p_[k] = slope(Vp_, Ip_);
    }

    this->setComponentBase(MEGA<RealT> * mva_);
  }

  template <class ScalarT, class IdxT>
  Reeca<ScalarT, IdxT>::~Reeca() = default;

  template <class ScalarT, class IdxT>
  int Reeca<ScalarT, IdxT>::allocate()
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
    if (auto& port = ports_.out[Output::iqcmd]; port)
      port.link(this->y_.getData() + IQCMD, &this->getVariableIndex(IQCMD));
    if (auto& port = ports_.out[Output::ipcmd]; port)
      port.link(this->y_.getData() + IPCMD, &this->getVariableIndex(IPCMD));
    this->allocated_ = true;
    return 0;
  }

  template <class ScalarT, class IdxT>
  GridKit::Model::ConfigurationChecks Reeca<ScalarT, IdxT>::verify() const
  {
    GridKit::Model::ConfigurationChecks checks("Reeca");
    checks.check(bus_ != nullptr, "terminal bus required");
    for (auto input : Utilities::enum_values<Input>())
    {
      if (ports_.in[input] && !ports_.in[input].linked())
        checks.fail(std::format("unlinked input {}", Utilities::enum_name(input)));
    }

    checks.check(mva_ > 0, "mva must be positive");
    checks.check(Trv_ >= 0 && Tp_ >= 0 && Tiq_ >= 0 && Tpord_ > 0, "invalid control time constants");
    checks.check(Vref0_ >= 0 && Vdip_ >= 0 && Vdip_ < Vup_, "invalid reference or dip voltages");
    checks.check(Thld2_ >= 0 && Imax_ > 0, "Thld2 must be nonnegative and Imax positive");
    checks.check(dbd1_ <= 0 && dbd2_ >= 0 && Iql1_ <= Iqh1_, "invalid injection limits or deadband");
    checks.check(dPmin_ < 0 && dPmax_ > 0 && Pmin_ <= Pmax_, "invalid power or ramp limits");
    checks.check(Vmin_ <= Vmax_ && Qmin_ <= Qmax_, "voltage and reactive limits must be ordered");
    const auto curve = [](const auto& v, const auto& i)
    {
      if (v == std::array<RealT, 4>{} && i == std::array<RealT, 4>{})
      {
        return true;
      }
      return v[0] >= 0 && v[0] < v[1] && v[1] < v[2] && v[2] < v[3] && std::ranges::min(i) >= 0;
    };
    checks.check(curve(Vq_, Iq_) && curve(Vp_, Ip_), "VDL tables must be empty or have increasing voltages and nonnegative limits");
    checks.check(this->ports_.out[Output::ipcmd] || this->ports_.out[Output::iqcmd], "current command output required");
    return checks;
  }

  template <class ScalarT, class IdxT>
  int Reeca<ScalarT, IdxT>::tagDifferentiable()
  {
    this->tag_[VMEAS] = Trv_ > 0;
    this->tag_[PMEAS] = Tp_ > 0;
    this->tag_[XPIQ]  = true;
    this->tag_[XPIV]  = true;
    this->tag_[QV]    = Tiq_ > 0;
    // The projected reference is algebraic while it tracks its input.
    this->tag_[PORD]  = true;
    return 0;
  }

  template <class ScalarT, class IdxT>
  template <class T>
  [[gnu::always_inline]] T Reeca<ScalarT, IdxT>::currentLimit(T voltage, const std::array<RealT, 4>& v, const std::array<RealT, 4>& i, const std::array<RealT, 3>& slope) const
  {
    T current{i.front() + (v.back() == 0) * Imax_};
#pragma clang loop unroll(full)
    for (size_t k = 0; k < 3; ++k)
    {
      current += slope[k] * (Math::ramp(voltage - v[k]) - Math::ramp(voltage - v[k + 1]));
    }
    return current - positiveDifference(current, Imax_);
  }

  template <class ScalarT, class IdxT>
  [[gnu::always_inline]] int Reeca<ScalarT, IdxT>::evaluateInternalResidual(const ScalarT* y, const ScalarT* yp, const ScalarT* wb, const ScalarT* w, ScalarT* f)
  {
    using T           = ScalarT;
    const auto signal = [w](Input input) -> const T&
    { return w[static_cast<size_t>(input)]; };

    const T speed    = ONE<RealT> + Pflag_ * (signal(Input::speed));
    const T pref     = this->toComponentBase(signal(Input::pref)) * speed;
    const T qref     = PfFlag_ * y[PMEAS] * std::tan(signal(Input::pfaref)) + !PfFlag_ * this->toComponentBase(signal(Input::qext));
    const T wanted   = Math::clamp(qref, Qmin_, Qmax_);
    const T vpi      = Kqp_ * y[EQ] + y[XPIQ];
    const T iqpi     = Kvp_ * y[EPIV] + y[XPIV];
    const T iq_table = currentLimit(y[VMEAS], Vq_, Iq_, slope_q_);
    const T ip_table = currentLimit(y[VMEAS], Vp_, Ip_, slope_p_);
    const T ipmax    = ip_table - !Pqflag_ * positiveDifference(ip_table, y[ILCAP]);
    const T iq       = QFlag_ * Math::clamp(iqpi, -y[IQMAX], y[IQMAX]) + !QFlag_ * y[QV]
                 + dip_ * y[IQV] + hold_q_ * ((Thld_ < 0) * y[IQV] + (Thld_ > 0) * Iqfrz_);
    const RealT run = ONE<RealT> - dip_;

    f[VSAFE] = Math::max(y[VMEAS], RealT{0.01}) - y[VSAFE];
    f[EQ]    = wanted - this->toComponentBase(signal(Input::qgen)) - y[EQ];
    f[EPIV]  = VFlag_ * Math::clamp(vpi, Vmin_, Vmax_)
              + !VFlag_ * (signal(Input::qext) + Vref1_) - y[VMEAS] - y[EPIV];
    f[ILCAP]     = GridKit::Model::currentCircle(Imax_, Pqflag_ * this->toComponentBase(y[IPCMD]) + !Pqflag_ * this->toComponentBase(y[IQCMD])) - y[ILCAP];
    f[IQMAX]     = iq_table - Pqflag_ * positiveDifference(iq_table, y[ILCAP]) - y[IQMAX];
    f[VTERM]     = std::sqrt(wb[0] * wb[0] + wb[1] * wb[1]) - y[VTERM];
    f[VMEAS]     = y[VTERM] - y[VMEAS] - Trv_ * yp[VMEAS];
    f[PMEAS]     = this->toComponentBase(signal(Input::pe)) - y[PMEAS] - Tp_ * yp[PMEAS];
    f[XPIQ]      = (QFlag_ && VFlag_) * run * Math::antiwindup(vpi, Kqi_ * y[EQ], Vmin_, Vmax_) - yp[XPIQ];
    f[XPIV]      = QFlag_ * run * Math::antiwindup(iqpi, Kvi_ * y[EPIV], -y[IQMAX], y[IQMAX]) - yp[XPIV];
    f[QV]        = (ONE<RealT> - (Tiq_ > 0) * dip_) * (wanted / y[VSAFE] - y[QV]) - Tiq_ * yp[QV];
    f[PREF_RATE] = std::clamp(pref - y[PREF_RATE] + yp[PREF_RATE], T{dPmin_}, T{dPmax_}) - yp[PREF_RATE];
    f[PORD]      = run * Math::antiwindup(y[PORD], (y[PREF_RATE] - y[PORD]) / Tpord_, Pmin_, Pmax_) - yp[PORD];
    f[IQV]       = Math::clamp(kqv_ * Math::deadband2(reference_voltage_ - y[VMEAS], dbd1_, dbd2_), Iql1_, Iqh1_) - y[IQV];
    f[IPMAX]     = (ONE<RealT> - hold_p_) * ipmax + hold_p_ * held_ipmax_ - y[IPMAX];
    f[IQCMD]     = Math::clamp(iq, -y[IQMAX], y[IQMAX]) - this->toComponentBase(y[IQCMD]);
    f[IPCMD]     = Math::clamp(y[PORD] / y[VSAFE], ZERO<RealT>, y[IPMAX]) - this->toComponentBase(y[IPCMD]);
    return 0;
  }

  template <class ScalarT, class IdxT>
  int Reeca<ScalarT, IdxT>::initialize()
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

    const RealT vt = std::hypot(static_cast<RealT>(bus_->Vr()), static_cast<RealT>(bus_->Vi()));
    if (vt <= Vdip_ || vt >= Vup_)
    {
      Log::error() << "Reeca: initial voltage outside the normal operating range\n";
      return 1;
    }
    const RealT safe = Math::max(vt, RealT{0.01});
    const RealT ip   = this->toComponentBase(ports_.out[Output::ipcmd] ? static_cast<RealT>(y[IPCMD]) : input(Input::pe) / safe);
    const RealT iq   = this->toComponentBase(ports_.out[Output::iqcmd] ? static_cast<RealT>(y[IQCMD]) : input(Input::qgen) / safe);
    y[IPCMD]         = this->toSystemBase(ip);
    y[IQCMD]         = this->toSystemBase(iq);
    pe_set_          = this->toSystemBase(ip * safe);
    qgen_set_        = this->toSystemBase(iq * safe);
    const RealT pe = this->toComponentBase(input(Input::pe, pe_set_)), q = this->toComponentBase(input(Input::qgen, qgen_set_));
    reference_voltage_ = Vref0_ > 0 ? Vref0_ : vt;
    dip_ = hold_p_ = hold_q_ = 0;
    const RealT circle       = GridKit::Model::currentCircle(Imax_, Pqflag_ ? ip : iq);
    const RealT iqmax        = std::min(currentLimit(vt, Vq_, Iq_, slope_q_), Pqflag_ ? circle : Imax_);
    const RealT ipmax        = std::min(currentLimit(vt, Vp_, Ip_, slope_p_), Pqflag_ ? Imax_ : circle);
    const auto  inverseClamp = [](RealT target, RealT lo, RealT hi)
    {
      return invertControl([&](RealT x)
                           { return Math::clamp(x, lo, hi); },
                           target,
                           lo - 1,
                           hi + 1);
    };
    const RealT pord    = safe * inverseClamp(ip, 0, ipmax);
    const RealT iqraw   = inverseClamp(iq, -iqmax, iqmax);
    const RealT command = QFlag_ ? q : safe * iqraw;
    const RealT qref    = inverseClamp(command, Qmin_, Qmax_);
    RealT       qext    = this->toSystemBase(qref);
    if (QFlag_ && !VFlag_)
      qext = vt - Vref1_;
    qext_set_ = qext;
    if (auto& port = ports_.in[Input::qext]; port)
      port.writeValue(ScalarT{qext_set_});
    pfaref_set_ = pe == 0 ? 0 : std::atan(qref / pe);
    if (auto& port = ports_.in[Input::pfaref]; port)
      port.writeValue(ScalarT{pfaref_set_});
    pref_set_ = this->toSystemBase(pord) / (1 + Pflag_ * input(Input::speed));
    if (auto& port = ports_.in[Input::pref]; port)
      port.writeValue(ScalarT{pref_set_});
    const RealT wanted = Math::clamp(PfFlag_ ? pe * std::tan(pfaref_set_) : this->toComponentBase(qext), Qmin_, Qmax_);
    const RealT eq     = wanted - q;
    const RealT vpi    = QFlag_ && VFlag_ ? inverseClamp(vt, Vmin_, Vmax_) : vt;
    const RealT epiv   = (VFlag_ ? Math::clamp(vpi, Vmin_, Vmax_) : qext + Vref1_) - vt;
    y[VTERM] = y[VMEAS] = vt;
    y[VSAFE]            = safe;
    y[PMEAS]            = pe;
    y[EQ]               = eq;
    y[EPIV]             = epiv;
    y[XPIQ]             = vpi - Kqp_ * eq;
    y[XPIV]             = (QFlag_ ? inverseClamp(iqraw, -iqmax, iqmax) : 0) - Kvp_ * epiv;
    y[QV]               = wanted / safe;
    y[PREF_RATE] = y[PORD] = pord;
    y[IQV]                 = Math::clamp(kqv_ * Math::deadband2(reference_voltage_ - vt, dbd1_, dbd2_), Iql1_, Iqh1_);
    y[ILCAP]               = circle;
    y[IQMAX]               = iqmax;
    y[IPMAX]               = ipmax;
    this->y_.setDataUpdated();
    evaluateResidual();
    for (size_t i = 0; i < COUNT; ++i)
    {
      const RealT r = static_cast<RealT>(this->f_.getData()[i]);
      if (!std::isfinite(r) || std::abs(r) > 1e-7 * (1 + std::abs(static_cast<RealT>(y[i]))))
      {
        Log::error() << "Reeca: infeasible initial operating point, row " << i << ", residual " << r << '\n';
        return 1;
      }
    }
    if constexpr (std::is_same_v<ScalarT, DependencyTracking::Variable>)
      this->initializeDependencyTrackingVariableNumbers();
    return 0;
  }

  template <class ScalarT, class IdxT>
  int Reeca<ScalarT, IdxT>::setSignalIndices()
  {
    for (auto input : Utilities::enum_values<Input>())
    {
      const auto& port                              = ports_.in[input];
      this->ws_indices_[static_cast<size_t>(input)] = port ? port.signalVariableIndex() : INVALID_INDEX<IdxT>;
    }
    return 0;
  }

  template <class ScalarT, class IdxT>
  int Reeca<ScalarT, IdxT>::evaluateResidual()
  {
    auto*      w      = this->ws_.getData();
    const auto gather = [&](Input input, RealT fallback = 0)
    {
      const auto& port              = ports_.in[input];
      w[static_cast<size_t>(input)] = port ? port.readSignal() : ScalarT{fallback};
    };
    gather(Input::pe, pe_set_);
    gather(Input::qgen, qgen_set_);
    gather(Input::qext, qext_set_);
    gather(Input::pfaref, pfaref_set_);
    gather(Input::pref, pref_set_);
    gather(Input::speed);
    this->wb_.getData()[0] = bus_->Vr();
    this->wb_.getData()[1] = bus_->Vi();
    evaluateInternalResidual(this->y_.getData(), this->yp_.getData(), this->wb_.getData(), w, this->f_.getData());
    this->f_.setDataUpdated();
    return 0;
  }

  template <typename ScalarT, typename IdxT>
  void Reeca<ScalarT, IdxT>::voltageDip(bool dip, RealT t, RealT ipmax)
  {
    if (dip == static_cast<bool>(dip_))
    {
      return;
    }
    dip_    = dip;
    hold_p_ = !dip && Thld2_ > 0;
    hold_q_ = !dip && Thld_ != 0;
    if (!dip)
    {
      held_ipmax_ = ipmax;
      release_p_  = t + Thld2_;
      release_q_  = t + std::abs(Thld_);
    }
  }

  template <typename ScalarT, typename IdxT>
  int Reeca<ScalarT, IdxT>::evaluateRootsAt(RealT t, const ScalarT* y, RealT* roots)
  {
    roots[0] = static_cast<RealT>(y[VTERM]) - Vdip_;
    roots[1] = static_cast<RealT>(y[VTERM]) - Vup_;
    roots[2] = hold_p_ != 0 ? t - release_p_ : RealT{1};
    roots[3] = hold_q_ != 0 ? t - release_q_ : RealT{1};
    return 0;
  }

  template <typename ScalarT, typename IdxT>
  bool Reeca<ScalarT, IdxT>::applyRootEvents(RealT t, std::span<const int> crossed)
  {
    const bool  was_dipped = dip_ != 0;
    const RealT ipmax      = static_cast<RealT>(this->y_.getData()[IPMAX]);
    if (crossed[2])
    {
      hold_p_ = 0;
    }
    if (crossed[3])
    {
      hold_q_ = 0;
    }
    if (crossed[0] || crossed[1])
    {
      voltageDip(crossed[0] < 0 || crossed[1] > 0, t, ipmax);
    }
    else
    {
      // Co-located controllers can reach the same threshold within IDA's root-time roundoff.
      const RealT v     = static_cast<RealT>(this->y_.getData()[VTERM]);
      const RealT dv    = static_cast<RealT>(this->yp_.getData()[VTERM]);
      const RealT round = 100 * std::numeric_limits<RealT>::epsilon() * (1 + std::abs(v) + std::abs(t * dv));
      if (dv != 0 && std::abs(v - Vdip_) <= round)
      {
        voltageDip(dv < 0, t, ipmax);
      }
      else if (dv != 0 && std::abs(v - Vup_) <= round)
      {
        voltageDip(dv > 0, t, ipmax);
      }
    }
    return was_dipped != (dip_ != 0) || crossed[2] || crossed[3];
  }

  template <class ScalarT, class IdxT>
  bool Reeca<ScalarT, IdxT>::synchronizeRootEvents(RealT t, const ScalarT* previous)
  {
    const RealT v   = static_cast<RealT>(this->y_.getData()[VTERM]);
    const bool  dip = v < Vdip_ || v > Vup_;
    if (dip == static_cast<bool>(dip_))
      return false;
    voltageDip(dip, t, static_cast<RealT>(previous[IPMAX]));
    return true;
  }
} // namespace GridKit::PhasorDynamics::Controller
