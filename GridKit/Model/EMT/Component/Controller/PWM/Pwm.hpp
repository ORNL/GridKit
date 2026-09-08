#pragma once

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Controller/PWM/PwmData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      /// Continuous pulse-width modulation with no DAE variables or residual rows.
      template <typename scalar_type, typename index_type>
      class Pwm : public Component<scalar_type, index_type>
      {
      public:
        using ScalarT    = scalar_type;
        using IdxT       = index_type;
        using RealT      = typename Component<ScalarT, IdxT>::RealT;
        using SignalT    = Signal<ScalarT, IdxT>;
        using ModelDataT = PwmData<RealT, IdxT>;
        using Inputs     = typename ModelDataT::Inputs;
        using Outputs    = typename ModelDataT::Outputs;
        using MonitorT   = Model::VariableMonitor<Pwm, PwmData>;

        Pwm();
        explicit Pwm(const ModelDataT& data);
        ~Pwm() override;

        int setGridKitComponentID(IdxT id) override final;
        int allocate() override final;
        int verify() const override final;

        int initialize(const std::map<Outputs, RealT>& outputs = {});

        typename Component<ScalarT, IdxT>::InitializationPortsT initializationPorts() override
        {
          return {};
        }

        int initializeState(const std::map<std::string, RealT>& values) override
        {
          return this->initializeOutputs(*this, values);
        }

        void validateInitialState(const std::map<std::string, RealT>& values) const override
        {
          this->template parseInitialOutputs<Pwm>(values);
        }

        int setAbsoluteTolerance(RealT) override final;
        int evaluateInternalResidual() override final;
        int evaluateExternalResidual() override final;
        int evaluateResidual() override final;
        int assembleJacobian(RealT y_scale, RealT yp_scale) override final;

        RealT maximumStepSize() const override final;

        void assignInput(Inputs key, SignalT* signal);

        /// Publish one phase on a named scalar signal. No DAE index is assigned.
        void    assignOutput(Outputs output, SignalT* signal);
        ScalarT output(Outputs output) const;

        SignalT& outputSignal(Outputs output)
        {
          return output_port_.at(static_cast<size_t>(output));
        }

      private:
        void                              initializeParameters(const ModelDataT& data);
        const Model::VariableMonitorBase* getMonitor() const override;
        bool                              hasInput() const;
        ScalarT                           modulation(size_t phase) const;
        ScalarT                           pulse(ScalarT duty, RealT local_time) const;
        void                              appendOutputGradient(Outputs output, typename SignalT::GradientT& gradient, RealT scale) const;

        RealT M_{0.0};
        RealT fm_{0.0};
        RealT fc_{0.0};
        RealT alignment_{0.5};
        bool  parameters_valid_{false};
        bool  sinusoidal_parameters_valid_{false};
        RealT horizon_{0.0};

        std::array<SignalT*, 3>   input_{};
        std::array<SignalT, 3>    output_port_;
        std::array<SignalT*, 3>   assigned_output_{};
        std::unique_ptr<MonitorT> monitor_;
      };
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
