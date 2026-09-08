#pragma once

#include <array>
#include <memory>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Controller/InnerCurrentControl/InnerCurrentControlData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      /// Synchronous-frame PI control with tracking anti-windup.
      template <typename scalar_type, typename index_type>
      class InnerCurrentControl : public Component<scalar_type, index_type>
      {
      public:
        using ScalarT    = scalar_type;
        using IdxT       = index_type;
        using Base       = Component<ScalarT, IdxT>;
        using RealT      = typename Base::RealT;
        using SignalT    = typename Base::SignalT;
        using ModelDataT = InnerCurrentControlData<RealT, IdxT>;
        using Inputs     = typename ModelDataT::Inputs;
        using Outputs    = typename ModelDataT::Outputs;
        using MonitorT   = Model::VariableMonitor<InnerCurrentControl, InnerCurrentControlData>;

        explicit InnerCurrentControl(const ModelDataT& data);
        ~InnerCurrentControl() override;

        int setGridKitComponentID(IdxT id) override final;
        int allocate() override final;
        int verify() const override final;

        int                                 initialize(const std::array<RealT, 2>& integral = {});
        void                                validateInitialState(const std::map<std::string, RealT>& values) const override;
        int                                 initializeState(const std::map<std::string, RealT>& values) override;
        typename Base::InitializationPortsT initializationPorts() override;

        int setAbsoluteTolerance(RealT tolerance) override final;
        int evaluateInternalResidual() override final;
        int evaluateResidual() override final;
        int assembleJacobian(RealT y_scale, RealT yp_scale) override final;

        void     attachInput(const std::array<SignalT*, 8>& inputs);
        void     assignOutput(Outputs output, SignalT* signal);
        SignalT& outputSignal(Outputs output);

      private:
        const Model::VariableMonitorBase* getMonitor() const override;
        ScalarT                           output(Outputs output) const;
        ScalarT                           dcVoltage() const;
        void                              appendOutputGradient(Outputs output, typename SignalT::GradientT& gradient, RealT scale) const;
        std::array<ScalarT, 2>            limitedReference() const;
        std::array<ScalarT, 2>            unlimitedVoltage() const;
        void                              appendLimitedGradient(size_t axis, typename SignalT::GradientT& gradient, RealT scale) const;
        void                              appendUnlimitedVoltageGradient(size_t axis, typename SignalT::GradientT& gradient, RealT scale) const;

        RealT                     inductance_, kp_, ki_, kaw_, current_limit_, modulation_limit_;
        RealT                     current_coefficient_, voltage_coefficient_;
        std::array<SignalT*, 8>   input_{};
        std::array<SignalT, 4>    output_{};
        std::array<SignalT*, 4>   alias_{};
        size_t                    capacity_{0};
        std::unique_ptr<MonitorT> monitor_;
      };
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
