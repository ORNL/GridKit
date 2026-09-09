/**
 * @file InnerCurrentControl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Inner-loop current controller.
 */
#pragma once

#include <array>
#include <memory>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Controller/InnerCurrentControl/InnerCurrentControlData.hpp>
#include <GridKit/Model/EMT/ComponentSignals.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      enum class InnerCurrentControlInternalVariables : size_t
      {
        XID,
        XIQ,
        ILIMD,
        ILIMQ,
        UD,
        UQ,
        MAXIMUM
      };
      enum class InnerCurrentControlExternalVariables : size_t
      {
        VD,
        VQ,
        ID,
        IQ,
        ICMDD,
        ICMDQ,
        OMEGA,
        ULIMD,
        ULIMQ,
        MAXIMUM
      };

      template <typename scalar_type, typename index_type>
      class InnerCurrentControl : public Component<scalar_type, index_type>
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
        using ModelDataT   = InnerCurrentControlData<RealT, IdxT>;
        using Outputs      = typename ModelDataT::Outputs;
        using SignalT      = Signal<ScalarT, IdxT>;
        using MonitorT     = Model::VariableMonitor<InnerCurrentControl, InnerCurrentControlData>;
        using InputSignals = std::array<SignalT*, 9>;

        InnerCurrentControl();
        explicit InnerCurrentControl(const ModelDataT& data);
        ~InnerCurrentControl();
        void     attachInput(InputSignals inputs);
        void     assignOutput(Outputs output, SignalT* signal);
        SignalT& outputSignal(Outputs output);

        SignalT& inputSignal(InnerCurrentControlInputs input)
        {
          return *signals_.getAttachedSignal(static_cast<InnerCurrentControlExternalVariables>(input));
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
        void                                                                                                        initializeParameters(const ModelDataT& data);
        void                                                                                                        initializeMonitor();
        const Model::VariableMonitorBase*                                                                           getMonitor() const override;
        RealT                                                                                                       L_{0.0};
        RealT                                                                                                       Kp_{0.0};
        RealT                                                                                                       Ki_{0.0};
        RealT                                                                                                       Kaw_{0.0};
        RealT                                                                                                       Imax_{0.0};
        RealT                                                                                                       ai_{0.0};
        RealT                                                                                                       i_scale_{ONE<RealT>};
        RealT                                                                                                       v_scale_{ONE<RealT>};
        ComponentSignals<ScalarT, IdxT, InnerCurrentControlInternalVariables, InnerCurrentControlExternalVariables> signals_;
        std::array<SignalT, 4>                                                                                      output_;
        std::array<SignalT*, 4>                                                                                     alias_{};
        std::unique_ptr<MonitorT>                                                                                   monitor_;
      };
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
