#pragma once

#include <cmath>

#include <GridKit/Model/EMT/Component/Controller/OuterVoltageControl/OuterVoltageControl.hpp>
#include <GridKit/Model/EMT/ComponentInitialization.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      template <typename scalar_type, typename index_type>
      OuterVoltageControl<scalar_type, index_type>::OuterVoltageControl()
        : OuterVoltageControl(ModelDataT{})
      {
      }

      template <typename scalar_type, typename index_type>
      OuterVoltageControl<scalar_type, index_type>::OuterVoltageControl(const ModelDataT& data)
        : monitor_(std::make_unique<MonitorT>(data))
      {
        initializeParameters(data);
        size_ = 4;
        signals_.template assignSignal<OuterVoltageControlInternalVariables::ICMDD>(&output_[0]);
        signals_.template assignSignal<OuterVoltageControlInternalVariables::ICMDQ>(&output_[1]);
        initializeMonitor();
      }

      template <typename scalar_type, typename index_type>
      OuterVoltageControl<scalar_type, index_type>::~OuterVoltageControl() = default;

      template <typename scalar_type, typename index_type>
      void OuterVoltageControl<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
      {
        using Parameter = typename ModelDataT::Parameters;
        C_              = parameter<RealT>(data, Parameter::C, C_);
        Kp_             = parameter<RealT>(data, Parameter::Kp, Kp_);
        Ki_             = parameter<RealT>(data, Parameter::Ki, Ki_);
        Kaw_            = parameter<RealT>(data, Parameter::Kaw, Kaw_);
      }

      template <typename scalar_type, typename index_type>
      void OuterVoltageControl<scalar_type, index_type>::attachInput(InputSignals inputs)
      {
        if (allocated_)
          throw std::logic_error("OuterVoltageControl inputs cannot change after allocation");
        for (size_t p = 0; p < inputs.size(); ++p)
          signals_.attachSignal(static_cast<OuterVoltageControlExternalVariables>(p), inputs[p]);
      }

      template <typename scalar_type, typename index_type>
      void OuterVoltageControl<scalar_type, index_type>::assignOutput(Outputs output, SignalT* signal)
      {
        if (allocated_)
          throw std::logic_error("OuterVoltageControl outputs cannot change after allocation");
        const auto n = static_cast<size_t>(output);
        if (n >= output_.size() || !signal || alias_[n])
          throw std::invalid_argument("OuterVoltageControl: invalid or duplicate output assignment");
        signal->claimProducer();
        alias_[n] = signal;
      }

      template <typename scalar_type, typename index_type>
      typename OuterVoltageControl<scalar_type, index_type>::SignalT& OuterVoltageControl<scalar_type, index_type>::outputSignal(Outputs output)
      {
        return output_.at(static_cast<size_t>(output));
      }

      template <typename scalar_type, typename index_type>
      int OuterVoltageControl<scalar_type, index_type>::setGridKitComponentID(IdxT id)
      {
        gridkit_component_id_ = id;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int OuterVoltageControl<scalar_type, index_type>::allocate()
      {
        if (!allocated_)
          this->allocateVectors(size_);
        tag_.resize(static_cast<size_t>(size_));
        variable_indices_.resize(static_cast<size_t>(size_));
        residual_indices_.resize(static_cast<size_t>(size_));
        this->assignGlobalIndices(0);
        this->allocateExternalVectors(static_cast<IdxT>(OuterVoltageControlExternalVariables::MAXIMUM), 0);
        signals_.registerExternalVariableSignals(*this);
        signals_.bindInternalVariableSignals(*this);
        for (size_t n = 0; n < alias_.size(); ++n)
          if (alias_[n])
            this->bindSignal(*alias_[n], static_cast<IdxT>(2 + n));
        allocated_ = true;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int OuterVoltageControl<scalar_type, index_type>::verify() const
      {
        int error_count   = 0;
        using V           = OuterVoltageControlExternalVariables;
        const auto inputs = signals_.attachedSignals({V::VREFD, V::VREFQ, V::VD, V::VQ, V::IGD, V::IGQ, V::OMEGA, V::ILIMD, V::ILIMQ});
        if (inputs.size() != static_cast<size_t>(V::MAXIMUM))
        {
          Log::error() << "OuterVoltageControl: all inputs are required\n";
          ++error_count;
        }
        for (const auto* input : inputs)
          if (!input->linked())
          {
            Log::error() << "OuterVoltageControl: inputs must have linked sources\n";
            ++error_count;
          }
        for (const auto value : {C_, Kp_, Ki_, Kaw_})
          if (!std::isfinite(value) || value <= ZERO<RealT>)
          {
            Log::error() << "OuterVoltageControl: parameters must be finite and positive\n";
            ++error_count;
          }
        return error_count;
      }

      template <typename scalar_type, typename index_type>
      void OuterVoltageControl<scalar_type, index_type>::validateInitialState(const std::map<std::string, RealT>& values) const
      {
        this->template parseInitialOutputs<OuterVoltageControl>(values);
      }

      template <typename scalar_type, typename index_type>
      int OuterVoltageControl<scalar_type, index_type>::initializeState(const std::map<std::string, RealT>& values)
      {
        return this->initializeOutputs(*this, values);
      }

      template <typename scalar_type, typename index_type>
      int OuterVoltageControl<scalar_type, index_type>::initialize(const std::map<Outputs, RealT>& outputs)
      {
        this->validateOutputValues(outputs);
        if (const int errors = verify(); errors != 0)
          return errors;
        using V           = OuterVoltageControlExternalVariables;
        const RealT vrefd = static_cast<RealT>(signals_.template readExternalVariable<V::VREFD>());
        const RealT vrefq = static_cast<RealT>(signals_.template readExternalVariable<V::VREFQ>());
        const RealT vd    = static_cast<RealT>(signals_.template readExternalVariable<V::VD>());
        const RealT vq    = static_cast<RealT>(signals_.template readExternalVariable<V::VQ>());
        const RealT igd   = static_cast<RealT>(signals_.template readExternalVariable<V::IGD>());
        const RealT igq   = static_cast<RealT>(signals_.template readExternalVariable<V::IGQ>());
        const RealT omega = static_cast<RealT>(signals_.template readExternalVariable<V::OMEGA>());
        for (const auto value : {vrefd, vrefq, vd, vq, igd, igq, omega})
          if (!std::isfinite(value))
            throw std::invalid_argument("OuterVoltageControl: nonfinite initial input");
        const RealT ed = vrefd - vd;
        const RealT eq = vrefq - vq;
        const RealT bd = igd - omega * C_ * vq;
        const RealT bq = igq + omega * C_ * vd;
        auto*       y  = y_.getData();
        y[2]           = this->outputValue(outputs, Outputs::icmdd, bd + Kp_ * ed);
        y[3]           = this->outputValue(outputs, Outputs::icmdq, bq + Kp_ * eq);
        y[0]           = y[2] - bd - Kp_ * ed;
        y[1]           = y[3] - bq - Kp_ * eq;
        auto* yp       = yp_.getData();
        for (IdxT n = 0; n < size_; ++n)
          yp[n] = ZERO<RealT>;
        y_.setDataUpdated();
        yp_.setDataUpdated();
        return 0;
      }

      template <typename scalar_type, typename index_type>
      typename Component<scalar_type, index_type>::InitializationPortsT OuterVoltageControl<scalar_type, index_type>::initializationPorts()
      {
        using V = OuterVoltageControlExternalVariables;
        // Tracking feedback is resolved by the consistent-initial-condition solve.
        return {signals_.attachedSignals({V::VREFD, V::VREFQ, V::VD, V::VQ, V::IGD, V::IGQ, V::OMEGA}), {}, {}};
      }

      template <typename scalar_type, typename index_type>
      int OuterVoltageControl<scalar_type, index_type>::setAbsoluteTolerance(RealT tolerance)
      {
        abs_tol_.setToConst(static_cast<ScalarT>(tolerance));
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int OuterVoltageControl<scalar_type, index_type>::evaluateInternalResidual(
          const ScalarT* y, const ScalarT* yp, const ScalarT* input, const ScalarT*, ScalarT* f)
      {
        const ScalarT ed = input[0] - input[2];
        const ScalarT eq = input[1] - input[3];
        const ScalarT bd = input[4] - input[6] * C_ * input[3];
        const ScalarT bq = input[5] + input[6] * C_ * input[2];
        f[0]             = -yp[0] + Ki_ * ed + Kaw_ * (input[7] - y[2]);
        f[1]             = -yp[1] + Ki_ * eq + Kaw_ * (input[8] - y[3]);
        f[2]             = y[2] - (bd + Kp_ * ed + y[0]);
        f[3]             = y[3] - (bq + Kp_ * eq + y[1]);
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int OuterVoltageControl<scalar_type, index_type>::evaluateInternalResidual()
      {
        this->gatherExternalVariables();
        evaluateInternalResidual(y_.getData(), yp_.getData(), y_ext_.data(), yp_ext_.data(), f_.getData());
        f_.setDataUpdated();
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int OuterVoltageControl<scalar_type, index_type>::evaluateResidual()
      {
        return evaluateInternalResidual();
      }

      template <typename scalar_type, typename index_type>
      void OuterVoltageControl<scalar_type, index_type>::initializeMonitor()
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
      const Model::VariableMonitorBase* OuterVoltageControl<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
