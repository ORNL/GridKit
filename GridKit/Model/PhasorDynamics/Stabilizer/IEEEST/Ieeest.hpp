/**
 * @file Ieeest.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the IEEEST Power System Stabilizer.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>
#include <GridKit/Model/PhasorDynamics/Stabilizer/IEEEST/IeeestData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename scalar_type, typename index_type>
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
        VSS, ///< Limited stabilizer signal (model output)
        MAXIMUM,
      };

      /// External variables of a `Ieeest`
      enum class IeeestExternalVariables : size_t
      {
        U, ///< Stabilizer input signal
        MAXIMUM,
      };

      template <typename scalar_type, typename index_type>
      class Ieeest : public Component<scalar_type, index_type>
      {
        using Component<scalar_type, index_type>::gridkit_component_id_;
        using Component<scalar_type, index_type>::alpha_;
        using Component<scalar_type, index_type>::f_;
        using Component<scalar_type, index_type>::nnz_;
        using Component<scalar_type, index_type>::size_;
        using Component<scalar_type, index_type>::tag_;
        using Component<scalar_type, index_type>::abs_tol_;
        using Component<scalar_type, index_type>::time_;
        using Component<scalar_type, index_type>::y_;
        using Component<scalar_type, index_type>::yp_;
        using Component<scalar_type, index_type>::wb_;
        using Component<scalar_type, index_type>::h_;
        using Component<scalar_type, index_type>::J_rows_buffer_;
        using Component<scalar_type, index_type>::J_cols_buffer_;
        using Component<scalar_type, index_type>::J_vals_buffer_;
        using Component<scalar_type, index_type>::variable_indices_;
        using Component<scalar_type, index_type>::residual_indices_;

      public:
        using ScalarT    = scalar_type;
        using IdxT       = index_type;
        using RealT      = typename Component<ScalarT, IdxT>::RealT;
        using ModelDataT = IeeestData<RealT, IdxT>;
        using SignalT    = SignalNode<ScalarT, IdxT>;
        using MonitorT   = Model::VariableMonitor<Ieeest, IeeestData>;

        Ieeest();
        Ieeest(const ModelDataT& data);
        ~Ieeest();

        int setGridKitComponentID(IdxT) override final;
        int allocate() override final;
        int verify() const override final;
        int initialize() override final;
        int tagDifferentiable() override final;
        int setAbsoluteTolerance(RealT rel_tol) override final;
        int evaluateResidual() override final;
        int evaluateJacobian() override final;

        /// Get the `ComponentSignals` from this `Ieeest`
        auto getSignals()
            -> ComponentSignals<ScalarT,
                                IdxT,
                                IeeestInputPorts,
                                IeeestOutputPorts>&
        {
          return signals_;
        }

        const Model::VariableMonitorBase* getMonitor() const override;

        __attribute__((always_inline)) inline int evaluateInternalResidual(
            const ScalarT*,
            const ScalarT*,
            const ScalarT*,
            const ScalarT*,
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
        RealT Vcu_{0};
        RealT Tdelay_{0};

        RealT a0_{1};
        RealT a1_{0};
        RealT a2_{0};
        RealT a3_{0};
        RealT a4_{0};

        // Precomputed masks and safe inverse coefficients for branch-free degenerate paths.
        RealT use_notch_{0};
        RealT bypass_notch_{1};
        RealT use_4th_order_{0};
        RealT use_3rd_order_{0};
        RealT use_2nd_order_{0};
        RealT safe_inv_a4_{0};
        RealT safe_inv_a3_{0};
        RealT safe_inv_a2_{0};
        RealT use_T2_block_{1};
        RealT bypass_T2_block_{0};
        RealT use_T4_block_{1};
        RealT bypass_T4_block_{0};
        RealT use_T6_block_{1};
        RealT bypass_T6_block_{0};

        ComponentSignals<ScalarT, IdxT, IeeestInputPorts, IeeestOutputPorts> signals_;

        std::unique_ptr<MonitorT> monitor_;

        void initializeParameters(const ModelDataT& data);
        void initializeMonitor();

        std::vector<ScalarT> ws_;
        std::vector<IdxT>    ws_indices_;
      };

    } // namespace Stabilizer
  } // namespace PhasorDynamics
} // namespace GridKit
