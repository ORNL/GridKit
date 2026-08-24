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
      template <typename real_type, typename index_type>
      struct IeeestData;
    } // namespace Stabilizer

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
      /// Internal variables of `Ieeest`
      enum class IeeestInternalVariables : size_t
      {
        X1,  ///< \f$x_1\f$ Notch-filter signal state
        X2,  ///< \f$x_2\f$ First derivative of the filtered signal
        X3,  ///< \f$x_3\f$ Second derivative of the filtered signal
        X4,  ///< \f$x_4\f$ Third derivative of the filtered signal
        X5,  ///< \f$x_5\f$ Lead-lag 1 state
        X6,  ///< \f$x_6\f$ Lead-lag 2 state
        X7,  ///< \f$x_7\f$ Washout state
        V4,  ///< \f$v_4\f$ Notch-filter output
        V5,  ///< \f$v_5\f$ Lead-lag 1 output
        V6,  ///< \f$v_6\f$ Lead-lag 2 output
        V7,  ///< \f$v_7\f$ Unlimited stabilizer signal
        VSS, ///< \f$V_{\mathrm{ss}}\f$ Limited stabilizer signal and model output
        MAXIMUM,
      };

      /// External variables of an `Ieeest`.
      enum class IeeestExternalVariables : size_t
      {
        U, ///< \f$u\f$ Stabilizer input signal
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
        using Component<scalar_type, index_type>::ws_;
        using Component<scalar_type, index_type>::ws_indices_;
        using Component<scalar_type, index_type>::h_;
        using Component<scalar_type, index_type>::J_rows_buffer_;
        using Component<scalar_type, index_type>::J_cols_buffer_;
        using Component<scalar_type, index_type>::J_vals_buffer_;
        using Component<scalar_type, index_type>::variable_indices_;
        using Component<scalar_type, index_type>::residual_indices_;
        using Component<scalar_type, index_type>::allocated_;

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
                                IeeestInternalVariables,
                                IeeestExternalVariables>&
        {
          return signals_;
        }

        const Model::VariableMonitorBase* getMonitor() const override;

        __attribute__((always_inline)) inline int evaluateInternalResidualKernel(
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

        ComponentSignals<ScalarT, IdxT, IeeestInternalVariables, IeeestExternalVariables> signals_;

        std::unique_ptr<MonitorT> monitor_;

        void initializeParameters(const ModelDataT& data);
        void initializeMonitor();
      };

    } // namespace Stabilizer
  } // namespace PhasorDynamics
} // namespace GridKit
