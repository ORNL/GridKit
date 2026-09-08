#pragma once

#include <array>
#include <memory>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Controller/OuterVoltageControl/OuterVoltageControlData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      /// Synchronous-frame PI control with tracking anti-windup.
      template <typename scalar_type, typename index_type>
      class OuterVoltageControl : public Component<scalar_type, index_type>
      {
      public:
        using ScalarT    = scalar_type;
        using IdxT       = index_type;
        using Base       = Component<ScalarT, IdxT>;
        using RealT      = typename Base::RealT;
        using SignalT    = typename Base::SignalT;
        using ModelDataT = OuterVoltageControlData<RealT, IdxT>;
        using Inputs     = typename ModelDataT::Inputs;
        using Outputs    = typename ModelDataT::Outputs;
        using MonitorT   = Model::VariableMonitor<OuterVoltageControl, OuterVoltageControlData>;

        explicit OuterVoltageControl(const ModelDataT& data);
        ~OuterVoltageControl() override;

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

        void     attachInput(const std::array<SignalT*, 9>& inputs);
        void     assignOutput(Outputs output, SignalT* signal);
        SignalT& outputSignal(Outputs output);

      private:
        const Model::VariableMonitorBase* getMonitor() const override;
        ScalarT                           output(Outputs output) const;
        void                              appendOutputGradient(Outputs output, typename SignalT::GradientT& gradient, RealT scale) const;

        RealT                     capacitance_, kp_, ki_, kaw_;
        std::array<SignalT*, 9>   input_{};
        std::array<SignalT, 2>    output_{};
        std::array<SignalT*, 2>   alias_{};
        size_t                    capacity_{0};
        std::unique_ptr<MonitorT> monitor_;
      };
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
