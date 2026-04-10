/**
 * @file Ieeest.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the IEEEST Power System Stabilizer.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Stabilizer
    {
      template <typename RealT, typename IdxT>
      struct IeeestData;
    } // namespace Stabilizer

    template <class ScalarT, typename IdxT>
    class SignalNode;

  } // namespace PhasorDynamics
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Stabilizer
    {
      /// Internal variables of a `Ieeest`
      enum class IeeestInternalVariables : size_t
      {
        X1,  ///< Notch filter state 1
        X2,  ///< Notch filter state 2
        X3,  ///< Notch filter state 3
        X4,  ///< Notch filter state 4
        X5,  ///< Lead-lag 1 state
        X6,  ///< Lead-lag 2 state
        X7,  ///< Washout state
        V4,  ///< Notch filter output
        V5,  ///< Lead-lag 1 output
        V6,  ///< Lead-lag 2 output
        V7,  ///< Unlimited stabilizer signal
        VSS, ///< Limited stabilizer signal
        VS,  ///< Stabilizer output
        MAXIMUM,
      };

      /// External variables of a `Ieeest`
      enum class IeeestExternalVariables : size_t
      {
        U,   ///< Stabilizer input signal
        VCT, ///< Cutout signal
        MAXIMUM,
      };

      template <class ScalarT, typename IdxT>
      class Ieeest : public Component<ScalarT, IdxT>
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
        using Component<ScalarT, IdxT>::h_;
        using Component<ScalarT, IdxT>::J_;
        using Component<ScalarT, IdxT>::J_rows_buffer_;
        using Component<ScalarT, IdxT>::J_cols_buffer_;
        using Component<ScalarT, IdxT>::J_vals_buffer_;
        using Component<ScalarT, IdxT>::variable_indices_;
        using Component<ScalarT, IdxT>::residual_indices_;

      public:
        using RealT           = typename Component<ScalarT, IdxT>::RealT;
        using model_data_type = IeeestData<RealT, IdxT>;
        using signal_type     = SignalNode<ScalarT, IdxT>;
        using MonitorT        = Model::VariableMonitor<Ieeest, IeeestData>;

        Ieeest();
        Ieeest(const model_data_type& data);
        ~Ieeest();

        int setGridKitComponentID(IdxT) override final;
        int allocate() override final;
        int verify() const override final;
        int initialize() override final;
        int tagDifferentiable() override final;
        int evaluateResidual() override final;
        int evaluateJacobian() override final;

        /// Get the `ComponentSignals` from this `Ieeest`
        auto getSignals()
            -> ComponentSignals<ScalarT,
                                IdxT,
                                IeeestInternalVariables,
                                IeeestExternalVariables>&
        {
          return signals_;
        }

        const Model::VariableMonitorBase* getMonitor() const override;

        __attribute__((always_inline)) inline int evaluateInternalResidual(
            ScalarT*,
            ScalarT*,
            [[maybe_unused]] ScalarT*,
            ScalarT*,
            ScalarT*);

      private:
        RealT A1_{0};
        RealT A2_{0};
        RealT A3_{0};
        RealT A4_{0};
        RealT A5_{0};
        RealT A6_{0};
        RealT T1_{0};
        RealT T2_{1};
        RealT T3_{0};
        RealT T4_{1};
        RealT T5_{0};
        RealT T6_{1};
        RealT Ks_{1};
        RealT Lsmin_{-0.1};
        RealT Lsmax_{0.1};
        RealT Vcl_{0};
        RealT Vcu_{1.5};
        RealT Tdelay_{0};

        RealT a0_{1};
        RealT a1_{0};
        RealT a2_{0};
        RealT a3_{0};
        RealT a4_{0};

        // Precomputed masks and safe inverse coefficients for branch-free degenerate paths.
        RealT use_4th_order_{0};
        RealT use_3rd_order_{0};
        RealT use_2nd_order_{0};
        RealT safe_inv_a4_{0};
        RealT safe_inv_a3_{0};
        RealT safe_inv_a2_{0};
        RealT use_T2_block_{1};
        RealT bypass_T2_block_{0};
        RealT safe_inv_T2_{1};
        RealT use_T4_block_{1};
        RealT bypass_T4_block_{0};
        RealT safe_inv_T4_{1};
        RealT use_T6_block_{1};
        RealT bypass_T6_block_{0};
        RealT safe_inv_T6_{1};
        RealT use_cutout_{1};
        RealT bypass_cutout_{0};

        ComponentSignals<ScalarT, IdxT, IeeestInternalVariables, IeeestExternalVariables> signals_;

        std::unique_ptr<MonitorT> monitor_;

        void initializeParameters(const model_data_type& data);
        void initializeMonitor();

        std::vector<ScalarT> ws_;
        std::vector<IdxT>    ws_indices_;
      };

    } // namespace Stabilizer
  } // namespace PhasorDynamics
} // namespace GridKit
