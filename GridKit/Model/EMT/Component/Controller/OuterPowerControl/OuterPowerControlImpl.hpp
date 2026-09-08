#pragma once

#include <cmath>

#include <GridKit/Model/EMT/Component/Controller/OuterPowerControl/OuterPowerControl.hpp>
#include <GridKit/Model/EMT/ComponentInitialization.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      template <typename scalar_type, typename index_type>
      OuterPowerControl<scalar_type, index_type>::OuterPowerControl()
        : OuterPowerControl(ModelDataT{})
      {
      }

      template <typename scalar_type, typename index_type>
      OuterPowerControl<scalar_type, index_type>::OuterPowerControl(const ModelDataT& data)
        : monitor_(std::make_unique<MonitorT>(data))
      {
        initializeParameters(data);
        size_ = 4;
        signals_.template assignSignal<OuterPowerControlInternalVariables::ICMDD>(&output_[0]);
        signals_.template assignSignal<OuterPowerControlInternalVariables::ICMDQ>(&output_[1]);
        initializeMonitor();
      }

      template <typename scalar_type, typename index_type>
      OuterPowerControl<scalar_type, index_type>::~OuterPowerControl() = default;

      template <typename scalar_type, typename index_type>
      void OuterPowerControl<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
      {
        using Parameter = typename ModelDataT::Parameters;
        V_              = parameter<RealT>(data, Parameter::V, V_);
        Pref_           = parameter<RealT>(data, Parameter::Pref, Pref_);
        Qref_           = parameter<RealT>(data, Parameter::Qref, Qref_);
        Kp_             = parameter<RealT>(data, Parameter::Kp, Kp_);
        Ki_             = parameter<RealT>(data, Parameter::Ki, Ki_);
        Kaw_            = parameter<RealT>(data, Parameter::Kaw, Kaw_);
        if (V_ > ZERO<RealT>)
        {
          irefd_ = Pref_ / V_;
          irefq_ = -Qref_ / V_;
        }
      }

      template <typename scalar_type, typename index_type>
      void OuterPowerControl<scalar_type, index_type>::attachInput(InputSignals inputs)
      {
        if (allocated_)
          throw std::logic_error("OuterPowerControl inputs cannot change after allocation");
        for (size_t p = 0; p < inputs.size(); ++p)
          signals_.attachSignal(static_cast<OuterPowerControlExternalVariables>(p), inputs[p]);
      }

      template <typename scalar_type, typename index_type>
      void OuterPowerControl<scalar_type, index_type>::assignOutput(Outputs output, SignalT* signal)
      {
        if (allocated_)
          throw std::logic_error("OuterPowerControl outputs cannot change after allocation");
        const auto n = static_cast<size_t>(output);
        if (n >= output_.size() || !signal || alias_[n])
          throw std::invalid_argument("OuterPowerControl: invalid or duplicate output assignment");
        signal->claimProducer();
        alias_[n] = signal;
      }

      template <typename scalar_type, typename index_type>
      typename OuterPowerControl<scalar_type, index_type>::SignalT& OuterPowerControl<scalar_type, index_type>::outputSignal(Outputs output)
      {
        return output_.at(static_cast<size_t>(output));
      }

      template <typename scalar_type, typename index_type>
      int OuterPowerControl<scalar_type, index_type>::setGridKitComponentID(IdxT id)
      {
        gridkit_component_id_ = id;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int OuterPowerControl<scalar_type, index_type>::allocate()
      {
        if (!allocated_)
          this->allocateVectors(size_);
        tag_.resize(static_cast<size_t>(size_));
        variable_indices_.resize(static_cast<size_t>(size_));
        residual_indices_.resize(static_cast<size_t>(size_));
        this->assignGlobalIndices(0);
        this->allocateExternalVectors(static_cast<IdxT>(OuterPowerControlExternalVariables::MAXIMUM), 0);
        signals_.registerExternalVariableSignals(*this);
        signals_.bindInternalVariableSignals(*this);
        for (size_t n = 0; n < alias_.size(); ++n)
          if (alias_[n])
            this->bindSignal(*alias_[n], static_cast<IdxT>(2 + n));
        allocated_ = true;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int OuterPowerControl<scalar_type, index_type>::verify() const
      {
        int error_count   = 0;
        using V           = OuterPowerControlExternalVariables;
        const auto inputs = signals_.attachedSignals({V::ID, V::IQ, V::ILIMD, V::ILIMQ});
        if (inputs.size() != static_cast<size_t>(V::MAXIMUM))
        {
          Log::error() << "OuterPowerControl: all inputs are required\n";
          ++error_count;
        }
        for (const auto* input : inputs)
          if (!input->linked())
          {
            Log::error() << "OuterPowerControl: inputs must have linked sources\n";
            ++error_count;
          }
        for (const auto value : {V_, Kp_, Ki_, Kaw_})
          if (!std::isfinite(value) || value <= ZERO<RealT>)
          {
            Log::error() << "OuterPowerControl: voltage rating and gains must be finite and positive\n";
            ++error_count;
          }
        if (!std::isfinite(Pref_) || !std::isfinite(Qref_)
            || !std::isfinite(irefd_) || !std::isfinite(irefq_))
        {
          Log::error() << "OuterPowerControl: setpoints and derived current references must be finite\n";
          ++error_count;
        }
        return error_count;
      }

      template <typename scalar_type, typename index_type>
      void OuterPowerControl<scalar_type, index_type>::validateInitialState(const std::map<std::string, RealT>& values) const
      {
        this->template parseInitialOutputs<OuterPowerControl>(values);
      }

      template <typename scalar_type, typename index_type>
      int OuterPowerControl<scalar_type, index_type>::initializeState(const std::map<std::string, RealT>& values)
      {
        return this->initializeOutputs(*this, values);
      }

      template <typename scalar_type, typename index_type>
      int OuterPowerControl<scalar_type, index_type>::initialize(const std::map<Outputs, RealT>& outputs)
      {
        this->validateOutputValues(outputs);
        if (const int errors = verify(); errors != 0)
          return errors;
        using V        = OuterPowerControlExternalVariables;
        const RealT id = static_cast<RealT>(signals_.template readExternalVariable<V::ID>());
        const RealT iq = static_cast<RealT>(signals_.template readExternalVariable<V::IQ>());
        for (const auto value : {id, iq})
          if (!std::isfinite(value))
            throw std::invalid_argument("OuterPowerControl: nonfinite initial current measurement");
        const RealT ed = irefd_ - id;
        const RealT eq = irefq_ - iq;
        auto*       y  = y_.getData();
        auto*       yp = yp_.getData();
        y[2]           = this->outputValue(outputs, Outputs::icmdd, Kp_ * ed);
        y[3]           = this->outputValue(outputs, Outputs::icmdq, Kp_ * eq);
        y[0]           = y[2] - Kp_ * ed;
        y[1]           = y[3] - Kp_ * eq;
        for (IdxT n = 0; n < size_; ++n)
          yp[n] = ZERO<RealT>;
        y_.setDataUpdated();
        yp_.setDataUpdated();
        return 0;
      }

      template <typename scalar_type, typename index_type>
      typename Component<scalar_type, index_type>::InitializationPortsT OuterPowerControl<scalar_type, index_type>::initializationPorts()
      {
        using V = OuterPowerControlExternalVariables;
        // The tracking input is resolved with the connected current limiter
        // by the consistent-initial-condition solve, after the integral is set.
        return {signals_.attachedSignals({V::ID, V::IQ}), {}, {}};
      }

      template <typename scalar_type, typename index_type>
      int OuterPowerControl<scalar_type, index_type>::setAbsoluteTolerance(RealT tolerance)
      {
        abs_tol_.setToConst(static_cast<ScalarT>(tolerance));
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int OuterPowerControl<scalar_type, index_type>::evaluateInternalResidual(
          const ScalarT* y, const ScalarT* yp, const ScalarT* input, const ScalarT*, ScalarT* f)
      {
        const ScalarT ed = irefd_ - input[0];
        const ScalarT eq = irefq_ - input[1];
        f[0]             = -yp[0] + Ki_ * ed + Kaw_ * (input[2] - y[2]);
        f[1]             = -yp[1] + Ki_ * eq + Kaw_ * (input[3] - y[3]);
        f[2]             = y[2] - (Kp_ * ed + y[0]);
        f[3]             = y[3] - (Kp_ * eq + y[1]);
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int OuterPowerControl<scalar_type, index_type>::evaluateInternalResidual()
      {
        this->gatherExternalVariables();
        evaluateInternalResidual(y_.getData(), yp_.getData(), y_ext_.data(), yp_ext_.data(), f_.getData());
        f_.setDataUpdated();
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int OuterPowerControl<scalar_type, index_type>::evaluateResidual()
      {
        return evaluateInternalResidual();
      }

      template <typename scalar_type, typename index_type>
      void OuterPowerControl<scalar_type, index_type>::initializeMonitor()
      {
        using Mon = typename ModelDataT::MonitorableVariables;
        monitor_->set(Mon::etad, [this]
                      { return y_.getData()[0]; });
        monitor_->set(Mon::etaq, [this]
                      { return y_.getData()[1]; });
        monitor_->set(Mon::icmdd, [this]
                      { return y_.getData()[2]; });
        monitor_->set(Mon::icmdq, [this]
                      { return y_.getData()[3]; });
      }

      template <typename scalar_type, typename index_type>
      const Model::VariableMonitorBase* OuterPowerControl<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
