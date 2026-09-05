#pragma once

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
      /// Sinusoidal pulse-width modulation with no internal variables or residual rows.
      template <typename scalar_type, typename index_type>
      class Pwm : public Component<scalar_type, index_type>
      {
      public:
        using ScalarT    = scalar_type;
        using IdxT       = index_type;
        using RealT      = typename Component<ScalarT, IdxT>::RealT;
        using SignalT    = Signal<ScalarT, IdxT>;
        using Port3T     = Port3<ScalarT, IdxT>;
        using ModelDataT = PwmData<RealT, IdxT>;
        using MonitorT   = Model::VariableMonitor<Pwm, PwmData>;

        Pwm();
        explicit Pwm(const ModelDataT& data);
        ~Pwm() override;

        int setGridKitComponentID(IdxT id) override final;
        int allocate() override final;
        int verify() const override final;
        int initialize() override final;
        int tagDifferentiable() override final;
        int setAbsoluteTolerance(RealT) override final;
        int evaluateInternalResidual() override final;
        int evaluateExternalResidual() override final;
        int evaluateResidual() override final;
        int evaluateJacobian() override final;

        /// Publish one phase on a named scalar signal. No DAE index is assigned.
        void    assignOutput(size_t phase, SignalT* signal);
        ScalarT output(size_t phase) const;

        Port3T& switchingPort()
        {
          return output_port_;
        }

      private:
        void                              initializeParameters(const ModelDataT& data);
        const Model::VariableMonitorBase* getMonitor() const override;

        RealT M_{0.0};
        RealT fm_{0.0};
        RealT fc_{0.0};
        RealT alignment_{0.5};
        bool  parameters_valid_{false};
        RealT horizon_{0.0};

        // Output workspace, independent of the DAE state.
        mutable std::array<RealT, 3> cached_time_{
            std::numeric_limits<RealT>::quiet_NaN(),
            std::numeric_limits<RealT>::quiet_NaN(),
            std::numeric_limits<RealT>::quiet_NaN()};
        mutable std::array<RealT, 3> cached_output_{};

        Port3T                    output_port_;
        std::array<SignalT*, 3>   assigned_output_{};
        std::unique_ptr<MonitorT> monitor_;
      };
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
