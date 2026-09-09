#pragma once

#include <cmath>
#include <numbers>

#include <GridKit/Model/EMT/ComponentInitialization.hpp>
#include <GridKit/Model/EMT/Operators/Reference/PLL/Pll.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    Pll<scalar_type, index_type>::Pll()
      : Pll(ModelDataT{})
    {
    }

    template <typename scalar_type, typename index_type>
    Pll<scalar_type, index_type>::Pll(const ModelDataT& data)
      : monitor_(std::make_unique<MonitorT>(data))
    {
      initializeParameters(data);
      size_ = 3;
      signals_.template assignSignal<PllInternalVariables::THETA>(&output_[0]);
      signals_.template assignSignal<PllInternalVariables::OMEGA>(&output_[1]);
      initializeMonitor();
    }

    template <typename scalar_type, typename index_type>
    Pll<scalar_type, index_type>::~Pll() = default;

    template <typename scalar_type, typename index_type>
    void Pll<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
    {
      using Parameter = typename ModelDataT::Parameters;
      V_              = parameter<RealT>(data, Parameter::V, V_);
      freq_           = parameter<RealT>(data, Parameter::f, freq_);
      Kp_             = parameter<RealT>(data, Parameter::Kp, Kp_);
      Ki_             = parameter<RealT>(data, Parameter::Ki, Ki_);
      omega0_         = TWO<RealT> * std::numbers::pi_v<RealT> * freq_;
      beta_scale_     = std::sqrt(THREE<RealT>) / TWO<RealT>;
      if (V_ > ZERO<RealT>)
        projection_scale_ = std::sqrt(TWO<RealT> / THREE<RealT>) / V_;
    }

    template <typename scalar_type, typename index_type>
    void Pll<scalar_type, index_type>::attachInput(PhaseSignals voltage)
    {
      if (allocated_)
        throw std::logic_error("PLL inputs cannot change after allocation");
      for (size_t p = 0; p < voltage.size(); ++p)
        signals_.attachSignal(static_cast<PllExternalVariables>(p), voltage[p]);
    }

    template <typename scalar_type, typename index_type>
    void Pll<scalar_type, index_type>::assignOutput(Outputs output, SignalT* signal)
    {
      if (allocated_)
        throw std::logic_error("PLL outputs cannot change after allocation");
      const auto n = static_cast<size_t>(output);
      if (n >= output_.size() || !signal || alias_[n])
        throw std::invalid_argument("PLL: invalid or duplicate output assignment");
      signal->claimProducer();
      alias_[n] = signal;
    }

    template <typename scalar_type, typename index_type>
    typename Pll<scalar_type, index_type>::SignalT& Pll<scalar_type, index_type>::outputSignal(Outputs output)
    {
      return output_.at(static_cast<size_t>(output));
    }

    template <typename scalar_type, typename index_type>
    int Pll<scalar_type, index_type>::setGridKitComponentID(IdxT id)
    {
      gridkit_component_id_ = id;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Pll<scalar_type, index_type>::allocate()
    {
      if (!allocated_)
        this->allocateVectors(size_);
      tag_.resize(static_cast<size_t>(size_));
      variable_indices_.resize(static_cast<size_t>(size_));
      residual_indices_.resize(static_cast<size_t>(size_));
      this->assignGlobalIndices(0);
      this->allocateExternalVectors(static_cast<IdxT>(PllExternalVariables::MAXIMUM), 0);
      signals_.registerExternalVariableSignals(*this);
      signals_.bindInternalVariableSignals(*this);
      for (size_t n = 0; n < alias_.size(); ++n)
        if (alias_[n])
          this->bindSignal(*alias_[n], static_cast<IdxT>(2 * n));
      allocated_ = true;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Pll<scalar_type, index_type>::verify() const
    {
      int error_count = 0;
      if (!signals_.template isAttached<PllExternalVariables::VA>()
          || !signals_.template isAttached<PllExternalVariables::VB>()
          || !signals_.template isAttached<PllExternalVariables::VC>())
      {
        Log::error() << "PLL: all bus voltage inputs are required\n";
        ++error_count;
      }
      else if (!signals_.template isLinked<PllExternalVariables::VA>()
               || !signals_.template isLinked<PllExternalVariables::VB>()
               || !signals_.template isLinked<PllExternalVariables::VC>())
      {
        Log::error() << "PLL: bus voltage inputs must have linked sources\n";
        ++error_count;
      }
      for (const auto value : {V_, freq_, Kp_, Ki_})
        if (!std::isfinite(value) || value <= ZERO<RealT>)
        {
          Log::error() << "PLL: ratings and gains must be finite and positive\n";
          ++error_count;
        }
      if (!std::isfinite(omega0_) || !std::isfinite(projection_scale_))
      {
        Log::error() << "PLL: nonfinite derived parameter\n";
        ++error_count;
      }
      return error_count;
    }

    template <typename scalar_type, typename index_type>
    void Pll<scalar_type, index_type>::validateInitialState(const std::map<std::string, RealT>& values) const
    {
      this->template parseInitialOutputs<Pll>(values);
    }

    template <typename scalar_type, typename index_type>
    int Pll<scalar_type, index_type>::initializeState(const std::map<std::string, RealT>& values)
    {
      return this->initializeOutputs(*this, values);
    }

    template <typename scalar_type, typename index_type>
    int Pll<scalar_type, index_type>::initialize(const std::map<Outputs, RealT>& outputs)
    {
      this->validateOutputValues(outputs);
      if (const int errors = verify(); errors != 0)
        return errors;
      this->gatherExternalVariables();
      for (const auto& value : y_ext_)
        if (!std::isfinite(static_cast<RealT>(value)))
          throw std::invalid_argument("PLL: nonfinite initial bus voltage");
      const RealT va    = static_cast<RealT>(y_ext_[0]);
      const RealT vb    = static_cast<RealT>(y_ext_[1]);
      const RealT vc    = static_cast<RealT>(y_ext_[2]);
      const RealT alpha = (TWO<RealT> * va - vb - vc) / THREE<RealT>;
      const RealT beta  = (vb - vc) / std::sqrt(THREE<RealT>);
      if (!outputs.contains(Outputs::theta) && std::hypot(alpha, beta) == ZERO<RealT>)
        throw std::invalid_argument("PLL: zero initial voltage requires theta");
      const RealT theta = this->outputValue(outputs, Outputs::theta, std::atan2(beta, alpha));
      const RealT omega = this->outputValue(outputs, Outputs::omega, omega0_);
      auto*       y     = y_.getData();
      auto*       yp    = yp_.getData();
      y[0]              = theta;
      const auto vq     = quadrature(y[0], y_ext_.data());
      y[1]              = (omega - omega0_ - Kp_ * vq) / Ki_;
      y[2]              = omega;
      yp[0]             = omega;
      yp[1]             = vq;
      yp[2]             = ZERO<RealT>;
      y_.setDataUpdated();
      yp_.setDataUpdated();
      return 0;
    }

    template <typename scalar_type, typename index_type>
    typename Component<scalar_type, index_type>::InitializationPortsT Pll<scalar_type, index_type>::initializationPorts()
    {
      using V = PllExternalVariables;
      typename Component<ScalarT, IdxT>::InitializationPortsT ports;
      ports.inputs = signals_.attachedSignals({V::VA, V::VB, V::VC});
      for (size_t n = 0; n < output_.size(); ++n)
      {
        const auto name = std::string(magic_enum::enum_name(static_cast<Outputs>(n)));
        ports.outputs.emplace(name, &output_[n]);
        if (alias_[n])
          ports.outputs.emplace(name, alias_[n]);
      }
      return ports;
    }

    template <typename scalar_type, typename index_type>
    void Pll<scalar_type, index_type>::prepareInitialization(typename Component<ScalarT, IdxT>::InitialStateT& initial)
    {
      if (initial.omega() == ZERO<RealT>)
        return;
      const RealT va    = initial.value(inputSignal(PllInputs::va));
      const RealT vb    = initial.value(inputSignal(PllInputs::vb));
      const RealT vc    = initial.value(inputSignal(PllInputs::vc));
      const RealT alpha = (TWO<RealT> * va - vb - vc) / THREE<RealT>;
      const RealT beta  = (vb - vc) / std::sqrt(THREE<RealT>);
      if (std::hypot(alpha, beta) == ZERO<RealT>)
        throw std::invalid_argument("PLL: balanced initialization requires nonzero voltage");
      initial.provide(output_[0], std::atan2(beta, alpha));
      initial.provide(output_[1], initial.omega());
    }

    template <typename scalar_type, typename index_type>
    int Pll<scalar_type, index_type>::setAbsoluteTolerance(RealT tolerance)
    {
      auto* absolute = abs_tol_.getData();
      absolute[0]    = tolerance;
      absolute[1]    = tolerance * omega0_ / Ki_;
      absolute[2]    = tolerance * omega0_;
      abs_tol_.setDataUpdated();
      return 0;
    }

    template <typename scalar_type, typename index_type>
    scalar_type Pll<scalar_type, index_type>::quadrature(ScalarT theta, const ScalarT* voltage) const
    {
      // Factor the three-phase projection through stationary alpha and beta.
      const ScalarT alpha = voltage[0] - HALF<RealT> * (voltage[1] + voltage[2]);
      const ScalarT beta  = beta_scale_ * (voltage[1] - voltage[2]);
      return projection_scale_ * (std::cos(theta) * beta - std::sin(theta) * alpha);
    }

    template <typename scalar_type, typename index_type>
    int Pll<scalar_type, index_type>::evaluateInternalResidual(
        const ScalarT* y, const ScalarT* yp, const ScalarT* voltage, const ScalarT*, ScalarT* f)
    {
      const ScalarT vq = quadrature(y[0], voltage);
      f[0]             = yp[0] - y[2];
      f[1]             = yp[1] - vq;
      f[2]             = y[2] - (omega0_ + Kp_ * vq + Ki_ * y[1]);
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Pll<scalar_type, index_type>::evaluateInternalResidual()
    {
      this->gatherExternalVariables();
      evaluateInternalResidual(y_.getData(), yp_.getData(), y_ext_.data(), yp_ext_.data(), f_.getData());
      f_.setDataUpdated();
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Pll<scalar_type, index_type>::evaluateResidual()
    {
      return evaluateInternalResidual();
    }

    template <typename scalar_type, typename index_type>
    void Pll<scalar_type, index_type>::initializeMonitor()
    {
      using Mon = typename ModelDataT::MonitorableVariables;
      monitor_->set(Mon::theta, [this]
                    { return y_.getData()[0]; });
      monitor_->set(Mon::xi, [this]
                    { return y_.getData()[1]; });
      monitor_->set(Mon::omega, [this]
                    { return y_.getData()[2]; });
      monitor_->set(Mon::vq, [this]
                    {
        this->gatherExternalVariables();
        return quadrature(y_.getData()[0], y_ext_.data()); });
    }

    template <typename scalar_type, typename index_type>
    const Model::VariableMonitorBase* Pll<scalar_type, index_type>::getMonitor() const
    {
      return monitor_.get();
    }
  } // namespace EMT
} // namespace GridKit
