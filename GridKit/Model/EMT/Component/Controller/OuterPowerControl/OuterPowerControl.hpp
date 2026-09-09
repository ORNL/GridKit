/**
 * @file OuterPowerControl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Outer power controller.
 */
#pragma once

#include <array>
#include <memory>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Controller/OuterPowerControl/OuterPowerControlData.hpp>
#include <GridKit/Model/EMT/ComponentSignals.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      enum class OuterPowerControlInternalVariables : size_t
      {
        ETAD,
        ETAQ,
        ICMDD,
        ICMDQ,
        MAXIMUM
      };
      enum class OuterPowerControlExternalVariables : size_t
      {
        VD,
        VQ,
        ID,
        IQ,
        ILIMD,
        ILIMQ,
        MAXIMUM
      };

      template <typename scalar_type, typename index_type>
      class OuterPowerControl : public Component<scalar_type, index_type>
      {
        using Component<scalar_type, index_type>::gridkit_component_id_;
        using Component<scalar_type, index_type>::size_;
        using Component<scalar_type, index_type>::nnz_;
        using Component<scalar_type, index_type>::time_;
        using Component<scalar_type, index_type>::alpha_;
        using Component<scalar_type, index_type>::y_;
        using Component<scalar_type, index_type>::yp_;
        using Component<scalar_type, index_type>::abs_tol_;
        using Component<scalar_type, index_type>::tag_;
        using Component<scalar_type, index_type>::y_ext_;
        using Component<scalar_type, index_type>::yp_ext_;
        using Component<scalar_type, index_type>::variable_indices_ext_;
        using Component<scalar_type, index_type>::f_;
        using Component<scalar_type, index_type>::J_rows_buffer_;
        using Component<scalar_type, index_type>::J_cols_buffer_;
        using Component<scalar_type, index_type>::J_vals_buffer_;
        using Component<scalar_type, index_type>::variable_indices_;
        using Component<scalar_type, index_type>::residual_indices_;
        using Component<scalar_type, index_type>::allocated_;

      public:
        using ScalarT      = scalar_type;
        using IdxT         = index_type;
        using RealT        = typename Component<ScalarT, IdxT>::RealT;
        using ModelDataT   = OuterPowerControlData<RealT, IdxT>;
        using Outputs      = typename ModelDataT::Outputs;
        using SignalT      = Signal<ScalarT, IdxT>;
        using MonitorT     = Model::VariableMonitor<OuterPowerControl, OuterPowerControlData>;
        using InputSignals = std::array<SignalT*, 6>;

        OuterPowerControl();
        explicit OuterPowerControl(const ModelDataT& data);
        ~OuterPowerControl();
        void     attachInput(InputSignals inputs);
        void     assignOutput(Outputs output, SignalT* signal);
        SignalT& outputSignal(Outputs output);

        SignalT& inputSignal(OuterPowerControlInputs input)
        {
          return *signals_.getAttachedSignal(static_cast<OuterPowerControlExternalVariables>(input));
        }

        int                                                     setGridKitComponentID(IdxT) override final;
        int                                                     allocate() override final;
        int                                                     verify() const override final;
        int                                                     initialize(const std::map<Outputs, RealT>& outputs = {});
        int                                                     initializeState(const std::map<std::string, RealT>& values) override;
        void                                                    validateInitialState(const std::map<std::string, RealT>& values) const override;
        typename Component<ScalarT, IdxT>::InitializationPortsT initializationPorts() override;
        void                                                    prepareInitialization(typename Component<ScalarT, IdxT>::InitialStateT& initial) override;
        int                                                     setAbsoluteTolerance(RealT) override final;
        int                                                     evaluateInternalResidual() override final;
        int                                                     evaluateResidual() override final;
        int                                                     assembleJacobian(RealT y_scale, RealT yp_scale) override final;

        auto& getSignals()
        {
          return signals_;
        }

        __attribute__((always_inline)) inline int evaluateInternalResidual(
            const ScalarT*, const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);

      private:
        void                                                                                                    initializeParameters(const ModelDataT& data);
        void                                                                                                    initializeMonitor();
        const Model::VariableMonitorBase*                                                                       getMonitor() const override;
        RealT                                                                                                   V_{0.0};
        RealT                                                                                                   Pref_{0.0};
        RealT                                                                                                   Qref_{0.0};
        RealT                                                                                                   Kp_{0.0};
        RealT                                                                                                   Ki_{0.0};
        RealT                                                                                                   Kaw_{0.0};
        ComponentSignals<ScalarT, IdxT, OuterPowerControlInternalVariables, OuterPowerControlExternalVariables> signals_;
        std::array<SignalT, 2>                                                                                  output_;
        std::array<SignalT*, 2>                                                                                 alias_{};
        std::unique_ptr<MonitorT>                                                                               monitor_;
      };
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
