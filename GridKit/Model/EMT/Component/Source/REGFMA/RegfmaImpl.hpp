#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/EMT/Component/Source/REGFMA/Regfma.hpp>
#include <GridKit/Model/EMT/ComponentInitialization.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    Regfma<scalar_type, index_type>::Regfma(const ModelDataT& data)
      : monitor_(std::make_unique<MonitorT>(data))
    {
      initializeParameters(data);
      this->equation_size_ = this->size_ = static_cast<IdxT>(I::MAXIMUM);
      for (size_t p = 0; p < 3; ++p)
        signals_.assignSignal(static_cast<I>(static_cast<size_t>(I::IA) + p), &current_[p]);
      initializeMonitor();
    }

    template <typename scalar_type, typename index_type>
    Regfma<scalar_type, index_type>::~Regfma() = default;

    template <typename scalar_type, typename index_type>
    void Regfma<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
    {
      using P          = typename ModelDataT::Parameters;
      S_               = parameter<RealT>(data, P::S);
      V_               = parameter<RealT>(data, P::V);
      omega0_          = parameter<RealT>(data, P::omega0, RealT{2} * std::numbers::pi_v<RealT> * RealT{60});
      XL_              = parameter<RealT>(data, P::XL, RealT{0.15});
      RL_              = parameter<RealT>(data, P::RL, RealT{0.03});
      mp_              = parameter<RealT>(data, P::mp, RealT{0.01});
      mq_              = parameter<RealT>(data, P::mq, RealT{0.05});
      kpv_             = parameter<RealT>(data, P::kpv, RealT{0});
      kiv_             = parameter<RealT>(data, P::kiv, RealT{5.86});
      Emin_            = parameter<RealT>(data, P::Emin, RealT{0});
      Emax_            = parameter<RealT>(data, P::Emax, RealT{1.15});
      Pmin_            = parameter<RealT>(data, P::Pmin, RealT{0});
      Pmax_            = parameter<RealT>(data, P::Pmax, RealT{0.9});
      Qmin_            = parameter<RealT>(data, P::Qmin, RealT{-0.44});
      Qmax_            = parameter<RealT>(data, P::Qmax, RealT{0.44});
      kppmax_          = parameter<RealT>(data, P::kppmax, RealT{0.01});
      kipmax_          = parameter<RealT>(data, P::kipmax, RealT{0.1});
      const bool VFlag = parameter<bool>(data, P::VFlag, true);
      kpqmax_          = parameter<RealT>(data, P::kpqmax, VFlag ? RealT{3} : RealT{0.1});
      kiqmax_          = parameter<RealT>(data, P::kiqmax, VFlag ? RealT{20} : RealT{10});
      TPf_             = parameter<RealT>(data, P::TPf, RealT{0.01});
      TQf_             = parameter<RealT>(data, P::TQf, RealT{0.01});
      TVf_             = parameter<RealT>(data, P::TVf, RealT{0.01});
      ImaxF_           = parameter<RealT>(data, P::ImaxF, RealT{2});
      QVFlag_          = parameter<bool>(data, P::QVFlag, true);
      voltage_control_ = VFlag ? RealT{1} : RealT{0};

      for (RealT value : {S_, V_, omega0_, XL_, RL_, mp_, TPf_, TQf_, TVf_, ImaxF_})
        if (!std::isfinite(value) || value <= RealT{0})
          throw std::invalid_argument("Regfma: bases, frequency, coupling impedance, active droop, filter times, and current limit must be positive");
      for (RealT value : {mq_, kpv_, kiv_, kppmax_, kipmax_, kpqmax_, kiqmax_})
        if (!std::isfinite(value) || value < RealT{0})
          throw std::invalid_argument("Regfma: control gains must be finite and nonnegative");
      for (RealT value : {Emin_, Emax_, Pmin_, Pmax_, Qmin_, Qmax_})
        if (!std::isfinite(value))
          throw std::invalid_argument("Regfma: limits must be finite");
      if (Emin_ < RealT{0} || Emin_ >= Emax_ || Pmin_ > Pmax_ || Qmin_ > Qmax_)
        throw std::invalid_argument("Regfma: invalid voltage or power limits");
      current_base_ = S_ / V_;
      kpP_          = kppmax_ / mp_;
      kiP_          = kipmax_ / mp_;
      if (!std::isfinite(current_base_) || current_base_ <= RealT{0}
          || !std::isfinite(kpP_) || !std::isfinite(kiP_)
          || !std::isfinite(RL_ * RL_ + XL_ * XL_) || RL_ * RL_ + XL_ * XL_ <= RealT{0}
          || !std::isfinite(omega0_ / XL_) || omega0_ / XL_ <= RealT{0})
        throw std::invalid_argument("Regfma: invalid derived base, gain, or coupling coefficient");
    }

    template <typename scalar_type, typename index_type>
    int Regfma<scalar_type, index_type>::setGridKitComponentID(IdxT id)
    {
      this->gridkit_component_id_ = id;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Regfma<scalar_type, index_type>::allocate()
    {
      if (!this->allocated_)
        this->allocateVectors(this->size_);
      this->tag_.resize(static_cast<size_t>(this->size_));
      this->variable_indices_.resize(static_cast<size_t>(this->size_));
      this->residual_indices_.resize(static_cast<size_t>(this->size_));
      this->assignGlobalIndices(0);
      this->allocateExternalVectors(static_cast<IdxT>(E::MAXIMUM), 0);
      signals_.registerExternalVariableSignals(*this);
      signals_.bindInternalVariableSignals(*this);
      for (size_t p = 0; p < 3; ++p)
        if (alias_[p])
          this->bindSignal(*alias_[p], static_cast<IdxT>(I::IA) + static_cast<IdxT>(p));
      this->allocated_ = true;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Regfma<scalar_type, index_type>::verify() const
    {
      if (!signals_.template isAttached<E::VA>() || !signals_.template isAttached<E::VB>()
          || !signals_.template isAttached<E::VC>())
        return 1;
      for (const auto* signal : signals_.attachedSignals({E::VA, E::VB, E::VC, E::PREF, E::QREF, E::VREF}))
        if (!signal->linked())
          return 1;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    void Regfma<scalar_type, index_type>::assignOutput(Outputs output, SignalT* signal)
    {
      if (!signal || output >= Outputs::SIZE)
        throw std::invalid_argument("Regfma: invalid current output");
      auto& assigned = alias_[static_cast<size_t>(output)];
      if (assigned == signal)
        return;
      if (assigned || this->allocated_)
        throw std::logic_error("Regfma: assign current outputs before allocation");
      signal->claimProducer();
      assigned = signal;
    }

    template <typename scalar_type, typename index_type>
    int Regfma<scalar_type, index_type>::initializeState(const std::map<std::string, RealT>& values)
    {
      return this->initializeOutputs(*this, values);
    }

    template <typename scalar_type, typename index_type>
    void Regfma<scalar_type, index_type>::validateInitialState(const std::map<std::string, RealT>& values) const
    {
      this->template parseInitialOutputs<Regfma>(values);
    }

    template <typename scalar_type, typename index_type>
    typename Regfma<scalar_type, index_type>::Base::InitializationPortsT
    Regfma<scalar_type, index_type>::initializationPorts()
    {
      typename Base::InitializationPortsT ports;
      ports.inputs = signals_.attachedSignals({E::VA, E::VB, E::VC, E::PREF, E::QREF, E::VREF});
      for (size_t p = 0; p < 3; ++p)
        ports.outputs.emplace(std::string("i") + "abc"[p], &current_[p]);
      return ports;
    }

    template <typename scalar_type, typename index_type>
    void Regfma<scalar_type, index_type>::gatherExternalVariables()
    {
      Base::gatherExternalVariables();
      if (!signals_.template isAttached<E::PREF>())
        this->y_ext_[static_cast<size_t>(E::PREF)] = pref_set_;
      if (!signals_.template isAttached<E::QREF>())
        this->y_ext_[static_cast<size_t>(E::QREF)] = qref_set_;
      if (!signals_.template isAttached<E::VREF>())
        this->y_ext_[static_cast<size_t>(E::VREF)] = vref_set_;
    }

    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline std::array<scalar_type, 6>
    Regfma<scalar_type, index_type>::measurements(const ScalarT* y, const ScalarT* ye) const
    {
      const RealT   a  = std::sqrt(RealT{2} / RealT{3});
      const RealT   b  = RealT{1} / std::numbers::sqrt2_v<RealT>;
      const ScalarT va = a * (ye[VA] - RealT{0.5} * (ye[VB] + ye[VC])) / V_;
      const ScalarT vb = b * (ye[VB] - ye[VC]) / V_;
      const ScalarT ia = a * (y[IA] - RealT{0.5} * (y[IB] + y[IC])) / current_base_;
      const ScalarT ib = b * (y[IB] - y[IC]) / current_base_;
      return {va, vb, ia, ib, va * ia + vb * ib, vb * ia - va * ib};
    }

    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline std::array<scalar_type, 5>
    Regfma<scalar_type, index_type>::controls(const ScalarT* y, const ScalarT* ye) const
    {
      const ScalarT p_limit = Math::min(kpP_ * (Pmax_ - y[PF]) + y[XPMAX], RealT{0})
                              + Math::max(kpP_ * (Pmin_ - y[PF]) + y[XPMIN], RealT{0});
      const ScalarT q_limit = Math::min(kpqmax_ * (Qmax_ - y[QF]) + y[XQMAX], RealT{0})
                              + Math::max(kpqmax_ * (Qmin_ - y[QF]) + y[XQMIN], RealT{0});
      const ScalarT command   = ye[VREF] + mq_ * (ye[QREF] - y[QF]) + q_limit;
      const ScalarT error     = command - y[VF];
      const ScalarT voltage   = Math::clamp((RealT{1} - voltage_control_) * command
                                              + voltage_control_ * (kpv_ * error + y[XV]),
                                          Emin_,
                                          Emax_);
      const ScalarT deviation = omega0_ * mp_ * (ye[PREF] - y[PF] + p_limit);
      return {p_limit, q_limit, error, voltage, deviation};
    }

    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline std::array<scalar_type, 2>
    Regfma<scalar_type, index_type>::source(const ScalarT* y, const ScalarT* ye) const
    {
      const RealT impedance_squared = RL_ * RL_ + XL_ * XL_;

      const auto    m       = measurements(y, ye);
      const auto    c       = controls(y, ye);
      const ScalarT theta   = omega0_ * this->time_ + y[DELTA];
      const ScalarT ea      = c[3] * std::cos(theta);
      const ScalarT eb      = c[3] * std::sin(theta);
      const ScalarT trial_a = (RL_ * (ea - m[0]) + XL_ * (eb - m[1])) / impedance_squared;
      const ScalarT trial_b = (RL_ * (eb - m[1]) - XL_ * (ea - m[0])) / impedance_squared;
      const ScalarT scale   = std::sqrt(Math::max(RealT{1}, (trial_a * trial_a + trial_b * trial_b) / (ImaxF_ * ImaxF_)));
      const ScalarT ia      = trial_a / scale;
      const ScalarT ib      = trial_b / scale;
      return {ia, ib};
    }

    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline int Regfma<scalar_type, index_type>::evaluateInternalResidual(
        const ScalarT* y, const ScalarT* yp, const ScalarT* ye, const ScalarT*, ScalarT* f)
    {
      const auto      m        = measurements(y, ye);
      const auto      c        = controls(y, ye);
      const auto      s        = source(y, ye);
      constexpr RealT infinity = std::numeric_limits<RealT>::infinity();
      f[PF]                    = -yp[PF] + (m[4] - y[PF]) / TPf_;
      f[QF]                    = -yp[QF] + (m[5] - y[QF]) / TQf_;
      f[VF]                    = -yp[VF] + (std::sqrt(m[0] * m[0] + m[1] * m[1] + VOLTAGE_EPSILON * VOLTAGE_EPSILON) - y[VF]) / TVf_;
      f[XPMAX]                 = -yp[XPMAX] + Math::antiwindup(y[XPMAX], kiP_ * (Pmax_ - y[PF]), -infinity, RealT{0});
      f[XPMIN]                 = -yp[XPMIN] + Math::antiwindup(y[XPMIN], kiP_ * (Pmin_ - y[PF]), RealT{0}, infinity);
      f[XQMAX]                 = -yp[XQMAX] + Math::antiwindup(y[XQMAX], kiqmax_ * (Qmax_ - y[QF]), -infinity, RealT{0});
      f[XQMIN]                 = -yp[XQMIN] + Math::antiwindup(y[XQMIN], kiqmax_ * (Qmin_ - y[QF]), RealT{0}, infinity);
      f[XV]                    = -yp[XV] + voltage_control_ * Math::antiwindup(y[XV], kiv_ * c[2], Emin_, Emax_);
      f[DELTA]                 = -yp[DELTA] + c[4];
      const RealT a            = std::sqrt(RealT{2} / RealT{3});
      const RealT b            = RealT{1} / std::numbers::sqrt2_v<RealT>;

      // L di/dt = e - v - R i, with the limited voltage drop (R_L + j X_L) i_lim.
      const ScalarT drop_a = current_base_ * (RL_ * s[0] - XL_ * s[1]);
      const ScalarT drop_b = current_base_ * (RL_ * s[1] + XL_ * s[0]);
      const RealT   rate   = omega0_ / XL_;
      f[IA]                = -yp[IA] + rate * (a * drop_a - RL_ * y[IA]);
      f[IB]                = -yp[IB] + rate * (-RealT{0.5} * a * drop_a + b * drop_b - RL_ * y[IB]);
      f[IC]                = -yp[IC] + rate * (-RealT{0.5} * a * drop_a - b * drop_b - RL_ * y[IC]);
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Regfma<scalar_type, index_type>::evaluateInternalResidual()
    {
      gatherExternalVariables();
      const int status = evaluateInternalResidual(this->y_.getData(), this->yp_.getData(), this->y_ext_.data(), this->yp_ext_.data(), this->f_.getData());
      this->f_.setDataUpdated();
      return status;
    }

    template <typename scalar_type, typename index_type>
    int Regfma<scalar_type, index_type>::evaluateResidual()
    {
      return evaluateInternalResidual();
    }

    template <typename scalar_type, typename index_type>
    typename Regfma<scalar_type, index_type>::RealT
    Regfma<scalar_type, index_type>::initialVoltageCommand(RealT voltage) const
    {
      if (!(voltage > Emin_ && voltage < Emax_))
        throw std::invalid_argument("Regfma: initial droop voltage must be strictly inside [Emin, Emax]");
      const RealT mu = Math::MU<RealT>;
      const RealT a  = mu * (voltage - Emin_);
      const RealT b  = mu * (Emax_ - voltage);
      return Emin_ + (a + std::log(-std::expm1(-a)) - std::log(-std::expm1(-b))) / mu;
    }

    template <typename scalar_type, typename index_type>
    int Regfma<scalar_type, index_type>::initialize(const std::map<Outputs, RealT>& outputs)
    {
      this->validateOutputValues(outputs);
      if (!this->allocated_ || verify() != 0)
        throw std::invalid_argument("Regfma: allocate and connect terminal voltages before initialization");
      gatherExternalVariables();
      for (const auto& value : this->y_ext_)
        if (!std::isfinite(static_cast<RealT>(value)))
          throw std::invalid_argument("Regfma: initial inputs must be finite");

      std::array<ScalarT, static_cast<size_t>(I::MAXIMUM)> state{};
      for (const auto& [output, value] : outputs)
        state[IA + static_cast<size_t>(output)] = value;
      const RealT sum = static_cast<RealT>(state[IA] + state[IB] + state[IC]);
      if (std::abs(sum) > RealT{1.0e-10} * current_base_)
        throw std::invalid_argument("Regfma: initial currents must have zero sum");

      const auto  m     = measurements(state.data(), this->y_ext_.data());
      const RealT va    = static_cast<RealT>(m[0]);
      const RealT vb    = static_cast<RealT>(m[1]);
      const RealT ia    = static_cast<RealT>(m[2]);
      const RealT ib    = static_cast<RealT>(m[3]);
      const RealT ratio = (ia * ia + ib * ib) / (ImaxF_ * ImaxF_);
      if (!(ratio < RealT{1}))
        throw std::invalid_argument("Regfma: initial current must be below ImaxF");

      // Invert the smooth radial limiter once, preserving the requested current.
      RealT lower = RealT{1};
      RealT upper = (RealT{1} + std::log(RealT{2}) / Math::MU<RealT>) / (RealT{1} - ratio);
      for (int iteration = 0; iteration < 100; ++iteration)
      {
        const RealT middle = lower + RealT{0.5} * (upper - lower);
        if (middle == lower || middle == upper)
          break;
        if (middle < Math::max(RealT{1}, ratio * middle))
          lower = middle;
        else
          upper = middle;
      }
      const RealT scale     = std::sqrt(lower + RealT{0.5} * (upper - lower));
      const RealT ea        = va + (RL_ * ia - XL_ * ib) * scale;
      const RealT eb        = vb + (RL_ * ib + XL_ * ia) * scale;
      const RealT magnitude = std::hypot(ea, eb);
      const RealT command   = initialVoltageCommand(magnitude);
      state[PF]             = m[4];
      state[QF]             = m[5];
      state[VF]             = std::sqrt(va * va + vb * vb + VOLTAGE_EPSILON * VOLTAGE_EPSILON);
      state[DELTA]          = std::atan2(eb, ea) - omega0_ * this->time_;

      auto input = this->y_ext_;
      if (!signals_.template isAttached<E::QREF>())
        input[QREF] = QVFlag_ ? ScalarT{0} : state[QF];
      const auto initial_controls = controls(state.data(), input.data());
      if (!signals_.template isAttached<E::PREF>())
        input[PREF] = state[PF] - initial_controls[0];
      if (!signals_.template isAttached<E::VREF>())
        input[VREF] = (voltage_control_ != RealT{0} ? state[VF] : ScalarT{command})
                      - mq_ * (input[QREF] - state[QF]) - initial_controls[1];
      const ScalarT error = input[VREF] + mq_ * (input[QREF] - state[QF]) + initial_controls[1] - state[VF];
      state[XV]           = voltage_control_ != RealT{0} ? ScalarT{command} - kpv_ * error : ScalarT{magnitude};
      if (voltage_control_ != RealT{0}
          && (static_cast<RealT>(state[XV]) < Emin_ || static_cast<RealT>(state[XV]) > Emax_))
        throw std::invalid_argument("Regfma: initial voltage integral is outside [Emin, Emax]");

      pref_set_ = input[PREF];
      qref_set_ = input[QREF];
      vref_set_ = input[VREF];
      std::copy(state.begin(), state.end(), this->y_.getData());
      this->yp_.setToConst(ScalarT{0});
      this->y_.setDataUpdated();
      evaluateResidual();
      for (size_t n = 0; n < static_cast<size_t>(I::MAXIMUM); ++n)
        this->yp_.getData()[n] = this->f_.getData()[n];
      this->yp_.setDataUpdated();
      return evaluateResidual();
    }

    template <typename scalar_type, typename index_type>
    int Regfma<scalar_type, index_type>::setAbsoluteTolerance(RealT tolerance)
    {
      this->abs_tol_.setToConst(static_cast<ScalarT>(tolerance));
      for (size_t p = 0; p < 3; ++p)
        this->abs_tol_.getData()[IA + p] = tolerance * current_base_;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    const Model::VariableMonitorBase* Regfma<scalar_type, index_type>::getMonitor() const
    {
      return monitor_.get();
    }

    template <typename scalar_type, typename index_type>
    void Regfma<scalar_type, index_type>::initializeMonitor()
    {
      using Mon = RegfmaMonitorableVariables;
      for (size_t n = 0; n < static_cast<size_t>(I::MAXIMUM); ++n)
        monitor_->set(static_cast<Mon>(n), [this, n]
                      { return this->y_.getData()[n]; });
      for (auto variable : {Mon::ea, Mon::eb, Mon::ec, Mon::omega, Mon::edroop, Mon::p, Mon::q, Mon::v})
        monitor_->set(variable, [this, variable]
                      { return monitorValue(variable); });
    }

    template <typename scalar_type, typename index_type>
    scalar_type Regfma<scalar_type, index_type>::monitorValue(RegfmaMonitorableVariables variable)
    {
      using Mon = RegfmaMonitorableVariables;
      gatherExternalVariables();
      const auto* y  = this->y_.getData();
      const auto* ye = this->y_ext_.data();
      if (variable == Mon::omega || variable == Mon::edroop)
      {
        const auto c = controls(y, ye);
        return variable == Mon::omega ? omega0_ + c[4] : c[3];
      }
      if (variable == Mon::p || variable == Mon::q || variable == Mon::v)
      {
        const auto m = measurements(y, ye);
        if (variable == Mon::v)
          return V_ * std::sqrt(m[0] * m[0] + m[1] * m[1] + VOLTAGE_EPSILON * VOLTAGE_EPSILON);
        return S_ * m[variable == Mon::p ? 4 : 5];
      }
      const auto    s  = source(y, ye);
      const ScalarT ea = RL_ * s[0] - XL_ * s[1];
      const ScalarT eb = RL_ * s[1] + XL_ * s[0];
      const RealT   a  = std::sqrt(RealT{2} / RealT{3});
      const RealT   b  = RealT{1} / std::numbers::sqrt2_v<RealT>;
      if (variable == Mon::ea)
        return ye[VA] + V_ * a * ea;
      return ye[variable == Mon::eb ? VB : VC] + V_ * (-RealT{0.5} * a * ea + (variable == Mon::eb ? b : -b) * eb);
    }
  } // namespace EMT
} // namespace GridKit
