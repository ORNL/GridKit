/**
 * @file SexsPti.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the SEXS-PTI exciter model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      template <typename RealT, typename IdxT>
      struct SexsPtiData;
    } // namespace Exciter

    template <class ScalarT, typename IdxT>
    class BusBase;

    template <class ScalarT, typename IdxT>
    class SignalNode;

  } // namespace PhasorDynamics
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      /// Internal variables of a `SexsPti`.
      enum class SexsPtiInternalVariables : size_t
      {
        VR,  ///< Lead-lag block state
        EFD, ///< Exciter field voltage output
        VTR, ///< Terminal voltage error signal
        MAXIMUM,
      };

      /// External variables of a `SexsPti`.
      enum class SexsPtiExternalVariables : size_t
      {
        VS, ///< Stabilizer output signal
        MAXIMUM,
      };

      template <class ScalarT, typename IdxT>
      class SexsPti : public Component<ScalarT, IdxT>
      {
        using Component<ScalarT, IdxT>::gridkit_component_id_;
        using Component<ScalarT, IdxT>::alpha_;
        using Component<ScalarT, IdxT>::f_;
        using Component<ScalarT, IdxT>::nnz_;
        using Component<ScalarT, IdxT>::size_;
        using Component<ScalarT, IdxT>::tag_;
        using Component<ScalarT, IdxT>::time_;
        using Component<ScalarT, IdxT>::y_;
        using Component<ScalarT, IdxT>::yp_;
        using Component<ScalarT, IdxT>::wb_;
        using Component<ScalarT, IdxT>::J_;
        using Component<ScalarT, IdxT>::J_rows_buffer_;
        using Component<ScalarT, IdxT>::J_cols_buffer_;
        using Component<ScalarT, IdxT>::J_vals_buffer_;
        using Component<ScalarT, IdxT>::variable_indices_;
        using Component<ScalarT, IdxT>::residual_indices_;

      public:
        using RealT           = typename Component<ScalarT, IdxT>::RealT;
        using model_data_type = SexsPtiData<RealT, IdxT>;
        using signal_type     = SignalNode<ScalarT, IdxT>;
        using bus_type        = BusBase<ScalarT, IdxT>;
        using MonitorT        = Model::VariableMonitor<SexsPti, SexsPtiData>;

        SexsPti(bus_type* bus);
        SexsPti(bus_type* bus, const model_data_type& data);
        ~SexsPti();

        int setGridKitComponentID(IdxT) override final;
        int allocate() override final;
        int verify() const override final;
        int initialize() override final;
        int tagDifferentiable() override final;
        int evaluateResidual() override final;
        int evaluateJacobian() override final;

        auto getSignals()
            -> ComponentSignals<ScalarT,
                                IdxT,
                                SexsPtiInternalVariables,
                                SexsPtiExternalVariables>&
        {
          return signals_;
        }

        const Model::VariableMonitorBase* getMonitor() const override;

        __attribute__((always_inline)) inline int evaluateInternalResidual(
            ScalarT*, ScalarT*, ScalarT*, ScalarT*, ScalarT*);

      private:
        bus_type* bus_{nullptr};

        RealT Ta_{0};
        RealT Tb_{0};
        RealT Te_{0};
        RealT K_{0};
        RealT Efdmax_{0};
        RealT Efdmin_{0};

        int missing_param_count_{0};

        ScalarT vref_{0};
        ScalarT vOEL_{0};
        ScalarT vUEL_{0};

        ComponentSignals<ScalarT, IdxT, SexsPtiInternalVariables, SexsPtiExternalVariables> signals_;

        std::unique_ptr<MonitorT> monitor_;

        void initModelParams(const model_data_type& data);
        void initializeMonitor();

        std::vector<ScalarT> ws_;
        std::vector<IdxT>    ws_indices_;
      };

    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
