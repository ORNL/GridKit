#pragma once

#include <deque>
#include <limits>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Controller/PWM/PwmData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      /// Sampled pulse-width modulation with no DAE variables or residual rows.
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
          // State owners initialize before this zero-state model samples its input.
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

        void  resetHistory() override final;
        void  acceptStep(RealT time) override final;
        RealT nextSampleTime(RealT after) const override final;
        RealT maximumStepSize() const override final;

        void assignInput(size_t phase, SignalT* signal);

        /// Publish one phase on a named scalar signal. No DAE index is assigned.
        void    assignOutput(size_t phase, SignalT* signal);
        ScalarT output(size_t phase) const;

        SignalT& outputSignal(Outputs output)
        {
          return output_port_.at(static_cast<size_t>(output));
        }

      private:
        void                              initializeParameters(const ModelDataT& data);
        const Model::VariableMonitorBase* getMonitor() const override;
        bool                              sampledInput() const;
        std::array<RealT, 3>              readModulation() const;
        void                              invalidateCache();
        RealT                             pulse(RealT duty, RealT local_time) const;
        RealT                             train(RealT duty, RealT t) const;
        RealT                             crossfade(RealT x) const;

        struct Sample
        {
          long long            interval;
          std::array<RealT, 3> modulation;
        };

        RealT              M_{0.0};
        RealT              fm_{0.0};
        RealT              fc_{0.0};
        RealT              alignment_{0.5};
        bool               parameters_valid_{false};
        bool               sinusoidal_parameters_valid_{false};
        RealT              horizon_{0.0};
        RealT              crossfade_rate_{0.0};
        std::deque<Sample> samples_;

        // Output workspace, independent of the DAE state.
        mutable std::array<RealT, 3> cached_time_{
            std::numeric_limits<RealT>::quiet_NaN(),
            std::numeric_limits<RealT>::quiet_NaN(),
            std::numeric_limits<RealT>::quiet_NaN()};
        mutable std::array<RealT, 3> cached_output_{};

        std::array<SignalT*, 3>   input_{};
        std::array<SignalT, 3>    output_port_;
        std::array<SignalT*, 3>   assigned_output_{};
        std::unique_ptr<MonitorT> monitor_;
      };
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
