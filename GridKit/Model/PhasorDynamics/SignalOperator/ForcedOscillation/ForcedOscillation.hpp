/**
 * @file ForcedOscillation.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the forced oscillation signal operator.
 */

#pragma once

#include <memory>
#include <vector>

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace SignalOperator
    {
      template <typename RealT, typename IdxT>
      struct ForcedOscillationData;
    } // namespace SignalOperator

    template <class ScalarT, typename IdxT>
    class SignalNode;
  } // namespace PhasorDynamics
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace SignalOperator
    {
      /// Internal signal variables owned by ForcedOscillation.
      enum class ForcedOscillationInternalVariables : size_t
      {
        OUT,    ///< Forced oscillation output signal
        MAXIMUM ///< Number of internal signal variables
      };

      /// External signal variables read by ForcedOscillation.
      enum class ForcedOscillationExternalVariables : size_t
      {
        IN,     ///< Optional input signal
        MAXIMUM ///< Number of external signal variables
      };

      /**
       * @brief Forced oscillation signal source/operator.
       *
       * @tparam ScalarT Scalar data type
       * @tparam IdxT    Matrix/vector index type
       */
      template <class ScalarT, typename IdxT>
      class ForcedOscillation : public Component<ScalarT, IdxT>
      {
        using Component<ScalarT, IdxT>::alpha_;
        using Component<ScalarT, IdxT>::f_;
        using Component<ScalarT, IdxT>::gridkit_component_id_;
        using Component<ScalarT, IdxT>::h_;
        using Component<ScalarT, IdxT>::J_;
        using Component<ScalarT, IdxT>::J_cols_buffer_;
        using Component<ScalarT, IdxT>::J_rows_buffer_;
        using Component<ScalarT, IdxT>::J_vals_buffer_;
        using Component<ScalarT, IdxT>::nnz_;
        using Component<ScalarT, IdxT>::residual_indices_;
        using Component<ScalarT, IdxT>::size_;
        using Component<ScalarT, IdxT>::tag_;
        using Component<ScalarT, IdxT>::time_;
        using Component<ScalarT, IdxT>::variable_indices_;
        using Component<ScalarT, IdxT>::wb_;
        using Component<ScalarT, IdxT>::y_;
        using Component<ScalarT, IdxT>::yp_;

      public:
        using RealT           = typename Component<ScalarT, IdxT>::RealT;
        using model_data_type = ForcedOscillationData<RealT, IdxT>;
        using signal_type     = SignalNode<ScalarT, IdxT>;
        using MonitorT        = Model::VariableMonitor<ForcedOscillation, ForcedOscillationData>;

        ForcedOscillation();
        ForcedOscillation(const model_data_type& data);
        ~ForcedOscillation();

        int  setGridKitComponentID(IdxT) override final;
        int  allocate() override final;
        int  verify() const override final;
        int  initialize() override final;
        void updateTime(RealT t, RealT a) override final;
        void initializeInputFromOutput();
        int  tagDifferentiable() override final;
        int  evaluateResidual() override final;
        int  evaluateJacobian() override final;

        /// Get the ComponentSignals owned by this ForcedOscillation.
        auto getSignals()
            -> ComponentSignals<ScalarT,
                                IdxT,
                                ForcedOscillationInternalVariables,
                                ForcedOscillationExternalVariables>&
        {
          return signals_;
        }

        const Model::VariableMonitorBase* getMonitor() const override;

        __attribute__((always_inline)) inline int evaluateInternalResidual(
            ScalarT*,
            [[maybe_unused]] ScalarT*,
            ScalarT*,
            ScalarT*,
            ScalarT*);

      private:
        void initializeParameters(const model_data_type& data);
        void initializeMonitor();
        void refreshForcing(RealT time);
        void refreshInitialOutput();

        RealT readRealParameter(const model_data_type&               data,
                                typename model_data_type::Parameters key,
                                RealT                                fallback);

        RealT A_{0.0};
        RealT f_param_{0.0};
        RealT Kf_{0.0};
        RealT Phi_{0.0};
        RealT Bias_{0.0};
        RealT Kin_{1.0};
        RealT u0_{0.0};
        RealT Ton_{0.0};
        RealT Toff_{-1.0};
        RealT Tr_{0.0};
        RealT Tf_{0.0};
        RealT Kd_{0.0};
        RealT Lmin_{-1.0e6};
        RealT Lmax_{1.0e6};

        ScalarT in_{0.0};
        ScalarT env_{0.0};
        ScalarT force_{0.0};
        ScalarT vraw_{0.0};
        ScalarT active_{0.0};

        size_t parameter_error_count_{0};

        ComponentSignals<ScalarT,
                         IdxT,
                         ForcedOscillationInternalVariables,
                         ForcedOscillationExternalVariables>
            signals_;

        std::unique_ptr<MonitorT> monitor_;

        std::vector<ScalarT> ws_;
        std::vector<IdxT>    ws_indices_;
      };

    } // namespace SignalOperator
  } // namespace PhasorDynamics
} // namespace GridKit
