#pragma once
#include <numeric>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/ParameterReader.hpp>
#include <GridKit/Model/PhasorDynamics/ControlLimits.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

#include "Ieeeg1.hpp"

namespace GridKit::PhasorDynamics::Governor
{

  template <class ScalarT, class IdxT>
  Ieeeg1<ScalarT, IdxT>::Ieeeg1(const DataT& data)
    : monitor_(std::make_unique<MonitorT>(data))
  {
    this->size_ = COUNT;
    monitor_->set(DataT::MonitorableVariables::pmech_hp, [this]
                  { return this->y_.getData()[PM]; });
    monitor_->set(DataT::MonitorableVariables::pmech_lp, [this]
                  { return this->y_.getData()[PM_LP]; });
    monitor_->set(DataT::MonitorableVariables::pgv, [this]
                  { return this->y_.getData()[VALVE]; });
    using Parameters = typename DataT::Parameters;
    const GridKit::Model::ParameterReader reader(data, "Ieeeg1");
    mva_            = reader.readReal(Parameters::mva).value_or(mva_);
    Trate_          = reader.readReal(Parameters::Trate).value_or(Trate_);
    K_              = reader.readReal(Parameters::K).value_or(K_);
    T1_             = reader.readReal(Parameters::T1).value_or(T1_);
    T2_             = reader.readReal(Parameters::T2).value_or(T2_);
    T3_             = reader.readReal(Parameters::T3).value_or(T3_);
    Uo_             = reader.readReal(Parameters::Uo).value_or(Uo_);
    Uc_             = reader.readReal(Parameters::Uc).value_or(Uc_);
    Pmax_           = reader.readReal(Parameters::Pmax).value_or(Pmax_);
    Pmin_           = reader.readReal(Parameters::Pmin).value_or(Pmin_);
    T4_             = reader.readReal(Parameters::T4).value_or(T4_);
    K1_             = reader.readReal(Parameters::K1).value_or(K1_);
    K2_             = reader.readReal(Parameters::K2).value_or(K2_);
    T5_             = reader.readReal(Parameters::T5).value_or(T5_);
    K3_             = reader.readReal(Parameters::K3).value_or(K3_);
    K4_             = reader.readReal(Parameters::K4).value_or(K4_);
    T6_             = reader.readReal(Parameters::T6).value_or(T6_);
    K5_             = reader.readReal(Parameters::K5).value_or(K5_);
    K6_             = reader.readReal(Parameters::K6).value_or(K6_);
    T7_             = reader.readReal(Parameters::T7).value_or(T7_);
    K7_             = reader.readReal(Parameters::K7).value_or(K7_);
    K8_             = reader.readReal(Parameters::K8).value_or(K8_);
    db1_            = reader.readReal(Parameters::db1).value_or(db1_);
    db2_            = reader.readReal(Parameters::db2).value_or(db2_);
    Eps_            = reader.readReal(Parameters::Eps).value_or(Eps_);
    Gv_[0]          = reader.readReal(Parameters::Gv1).value_or(0);
    Gv_[1]          = reader.readReal(Parameters::Gv2).value_or(0);
    Gv_[2]          = reader.readReal(Parameters::Gv3).value_or(0);
    Gv_[3]          = reader.readReal(Parameters::Gv4).value_or(0);
    Gv_[4]          = reader.readReal(Parameters::Gv5).value_or(0);
    Gv_[5]          = reader.readReal(Parameters::Gv6).value_or(0);
    Pgv_[0]         = reader.readReal(Parameters::Pgv1).value_or(0);
    Pgv_[1]         = reader.readReal(Parameters::Pgv2).value_or(0);
    Pgv_[2]         = reader.readReal(Parameters::Pgv3).value_or(0);
    Pgv_[3]         = reader.readReal(Parameters::Pgv4).value_or(0);
    Pgv_[4]         = reader.readReal(Parameters::Pgv5).value_or(0);
    Pgv_[5]         = reader.readReal(Parameters::Pgv6).value_or(0);
    lead_gain_      = T1_ > 0 ? T2_ / T1_ : ZERO<RealT>;
    straight_valve_ = Gv_ == std::array<RealT, 6>{} && Pgv_ == std::array<RealT, 6>{};
    for (size_t k = 0; k < slope_.size(); ++k)
    {
      if (Gv_[k + 1] > Gv_[k])
      {
        slope_[k] = (Pgv_[k + 1] - Pgv_[k]) / (Gv_[k + 1] - Gv_[k]);
      }
    }

    this->setComponentBase(MEGA<RealT> * (Trate_ > 0 ? Trate_ : mva_));
  }

  template <class ScalarT, class IdxT>
  Ieeeg1<ScalarT, IdxT>::~Ieeeg1() = default;

  template <class ScalarT, class IdxT>
  int Ieeeg1<ScalarT, IdxT>::allocate()
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
    if (auto& port = ports_.out[Output::pmech_hp]; port)
      port.link(this->y_.getData() + PM, &this->getVariableIndex(PM));
    if (auto& port = ports_.out[Output::pmech_lp]; port)
      port.link(this->y_.getData() + PM_LP, &this->getVariableIndex(PM_LP));
    this->allocated_ = true;
    return 0;
  }

  template <class ScalarT, class IdxT>
  GridKit::Model::ConfigurationChecks Ieeeg1<ScalarT, IdxT>::verify() const
  {
    GridKit::Model::ConfigurationChecks checks("Ieeeg1");
    for (auto input : Utilities::enum_values<Input>())
    {
      if (ports_.in[input] && !ports_.in[input].linked())
        checks.fail(std::format("unlinked input {}", Utilities::enum_name(input)));
    }

    checks.check((Trate_ > 0 || mva_ > 0) && Trate_ >= 0, "Trate must be nonnegative and requires positive mva when zero");
    checks.check(T1_ >= 0 && T2_ >= 0 && (T1_ > 0 || T2_ == 0), "T2 requires a positive T1");
    checks.check(T3_ > 0 && T4_ >= 0 && T5_ >= 0 && T6_ >= 0 && T7_ >= 0, "invalid servo or stage time constant");
    checks.check(Uc_ < 0 && Uo_ > 0 && Pmin_ <= Pmax_, "invalid valve limits");
    checks.check(db1_ >= 0 && Eps_ == 0 && db2_ == 0, "db1 must be nonnegative; hysteresis and backlash require zero Eps and db2");
    if (Gv_ != decltype(Gv_){} || Pgv_ != decltype(Pgv_){})
    {
      for (size_t k = 1; k < Gv_.size(); ++k)
      {
        checks.check(Gv_[k] > Gv_[k - 1] && Pgv_[k] >= Pgv_[k - 1], "valve curve must be monotone with rising Gv");
      }
    }
    checks.check(this->ports_.out[Output::pmech_hp] || this->ports_.out[Output::pmech_lp], "a shaft output is required");
    return checks;
  }

  template <class ScalarT, class IdxT>
  int Ieeeg1<ScalarT, IdxT>::tagDifferentiable()
  {
    this->tag_[SPEED]     = T1_ > 0;
    this->tag_[VALVE]     = true;
    this->tag_[BOWL]      = T4_ > 0;
    this->tag_[REHEAT]    = T5_ > 0;
    this->tag_[CROSSOVER] = T6_ > 0;
    this->tag_[REHEAT2]   = T7_ > 0;
    return 0;
  }

  template <class ScalarT, class IdxT>
  template <class T>
  [[gnu::always_inline]] T Ieeeg1<ScalarT, IdxT>::valveFlow(T valve) const
  {
    T flow = straight_valve_ * valve + Pgv_.front();
#pragma clang loop unroll(full)
    for (size_t k = 0; k < slope_.size(); ++k)
    {
      flow += slope_[k] * (Math::ramp(valve - Gv_[k]) - Math::ramp(valve - Gv_[k + 1]));
    }

    return flow;
  }

  template <class ScalarT, class IdxT>
  [[gnu::always_inline]] int Ieeeg1<ScalarT, IdxT>::evaluateInternalResidual(const ScalarT* y, const ScalarT* yp, const ScalarT*, const ScalarT* w, ScalarT* f)
  {
    using T           = ScalarT;
    const auto signal = [w](Input input) -> const T&
    { return w[static_cast<size_t>(input)]; };

    const T deviation = signal(Input::speed);
    const T speed     = K_ * Math::deadband1(deviation, -db1_, db1_);
    const T lead      = y[SPEED] + lead_gain_ * (speed - y[SPEED]);
    const T rate      = (this->toComponentBase(signal(Input::pref) + signal(Input::paux)) - lead - y[VALVE]) / T3_;
    const T flow      = valveFlow(y[VALVE]);
    f[SPEED]          = speed - y[SPEED] - T1_ * yp[SPEED];
    f[VALVE]          = Math::boundedAntiwindup(y[VALVE], Math::clamp(rate, Uc_, Uo_), Pmin_, Pmax_) - yp[VALVE];
    f[BOWL]           = flow - y[BOWL] - T4_ * yp[BOWL];
    f[REHEAT]         = y[BOWL] - y[REHEAT] - T5_ * yp[REHEAT];
    f[CROSSOVER]      = y[REHEAT] - y[CROSSOVER] - T6_ * yp[CROSSOVER];
    f[REHEAT2]        = y[CROSSOVER] - y[REHEAT2] - T7_ * yp[REHEAT2];
    f[PM]             = K1_ * y[BOWL] + K3_ * y[REHEAT] + K5_ * y[CROSSOVER] + K7_ * y[REHEAT2] - this->toComponentBase(y[PM]);
    f[PM_LP]          = K2_ * y[BOWL] + K4_ * y[REHEAT] + K6_ * y[CROSSOVER] + K8_ * y[REHEAT2] - this->toComponentBase(y[PM_LP]);
    return 0;
  }

  template <class ScalarT, class IdxT>
  int Ieeeg1<ScalarT, IdxT>::initialize()
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

    const RealT kh = K1_ + K3_ + K5_ + K7_, kl = K2_ + K4_ + K6_ + K8_;
    const bool  use_hp = ports_.out[Output::pmech_hp] && kh != 0;
    if (!use_hp && (!ports_.out[Output::pmech_lp] || kl == 0))
      return 1;
    const RealT flow  = this->toComponentBase(use_hp ? static_cast<RealT>(y[PM]) / kh : static_cast<RealT>(y[PM_LP]) / kl);
    const RealT valve = invertControl([&](RealT x)
                                      { return valveFlow(x); },
                                      flow,
                                      Pmin_,
                                      Pmax_);
    y[BOWL] = y[REHEAT] = y[CROSSOVER] = y[REHEAT2] = flow;
    y[VALVE]                                        = valve;
    const RealT speed                               = K_ * Math::deadband1(input(Input::speed), -db1_, db1_);
    y[SPEED]                                        = speed;
    if (!ports_.out[Output::pmech_lp])
      y[PM_LP] = this->toSystemBase(kl * flow);
    if (!ports_.out[Output::pmech_hp])
      y[PM] = this->toSystemBase(kh * flow);
    const RealT steady_rate = invertControl([&](RealT x)
                                            { return Math::clamp(x, Uc_, Uo_); },
                                            RealT{0},
                                            Uc_ - 1,
                                            Uo_ + 1);
    pref_set_               = this->toSystemBase(valve + speed + T3_ * steady_rate) - input(Input::paux);
    if (auto& port = ports_.in[Input::pref]; port)
      port.writeValue(ScalarT{pref_set_});
    this->y_.setDataUpdated();
    evaluateResidual();
    for (size_t i = 0; i < COUNT; ++i)
    {
      const RealT r = static_cast<RealT>(this->f_.getData()[i]);
      if (!std::isfinite(r) || std::abs(r) > 1e-7 * (1 + std::abs(static_cast<RealT>(y[i]))))
      {
        Log::error() << "Ieeeg1: infeasible initial operating point, row " << i << ", residual " << r << '\n';
        return 1;
      }
    }
    if constexpr (std::is_same_v<ScalarT, DependencyTracking::Variable>)
      this->initializeDependencyTrackingVariableNumbers();
    return 0;
  }

  template <class ScalarT, class IdxT>
  int Ieeeg1<ScalarT, IdxT>::setSignalIndices()
  {
    for (auto input : Utilities::enum_values<Input>())
    {
      const auto& port                              = ports_.in[input];
      this->ws_indices_[static_cast<size_t>(input)] = port ? port.signalVariableIndex() : INVALID_INDEX<IdxT>;
    }
    return 0;
  }

  template <class ScalarT, class IdxT>
  int Ieeeg1<ScalarT, IdxT>::evaluateResidual()
  {
    auto*      w      = this->ws_.getData();
    const auto gather = [&](Input input, RealT fallback = 0)
    {
      const auto& port              = ports_.in[input];
      w[static_cast<size_t>(input)] = port ? port.readSignal() : ScalarT{fallback};
    };
    gather(Input::speed);
    gather(Input::pref, pref_set_);
    gather(Input::paux);
    evaluateInternalResidual(this->y_.getData(), this->yp_.getData(), nullptr, w, this->f_.getData());
    this->f_.setDataUpdated();
    return 0;
  }

} // namespace GridKit::PhasorDynamics::Governor
