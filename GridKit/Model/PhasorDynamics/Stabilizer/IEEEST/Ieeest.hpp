/**
 * @file Ieeest.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the IEEEST Power System Stabilizer.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ConnectedElement.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Stabilizer
    {
      template <typename RealP, typename IdxP>
      struct IeeestData;
    } // namespace Stabilizer
  } // namespace PhasorDynamics
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Stabilizer
    {
      template <typename, typename>
      class Ieeest;

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
    } // namespace Stabilizer

    template <typename ScalarP, typename IdxP>
    struct ConnectedElementTraits<Stabilizer::Ieeest<ScalarP, IdxP>>
    {
      using IeeestT = Stabilizer::Ieeest<ScalarP, IdxP>;

      using ElementT           = IeeestT;
      using ScalarT            = ScalarP;
      using IdxT               = IdxP;
      using RealT              = typename ScalarTraits<ScalarT>::RealT;
      using ModelDataT         = Stabilizer::IeeestData<RealT, IdxT>;
      using InternalVariablesT = Stabilizer::IeeestInternalVariables;
      using ExternalVariablesT = Stabilizer::IeeestExternalVariables;
      using InterfaceT         = Component<ScalarT, IdxT>;
    };

    namespace Stabilizer
    {
      template <typename ScalarP, typename IdxP>
      class Ieeest : public ConnectedElement<Ieeest<ScalarP, IdxP>>
      {
        using ConnectedElement<Ieeest>::gridkit_component_id_;
        using ConnectedElement<Ieeest>::alpha_;
        using ConnectedElement<Ieeest>::f_;
        using ConnectedElement<Ieeest>::nnz_;
        using ConnectedElement<Ieeest>::size_;
        using ConnectedElement<Ieeest>::tag_;
        using ConnectedElement<Ieeest>::time_;
        using ConnectedElement<Ieeest>::y_;
        using ConnectedElement<Ieeest>::yp_;
        using ConnectedElement<Ieeest>::wb_;
        using ConnectedElement<Ieeest>::h_;
        using ConnectedElement<Ieeest>::J_;
        using ConnectedElement<Ieeest>::J_rows_buffer_;
        using ConnectedElement<Ieeest>::J_cols_buffer_;
        using ConnectedElement<Ieeest>::J_vals_buffer_;
        using ConnectedElement<Ieeest>::variable_indices_;
        using ConnectedElement<Ieeest>::residual_indices_;
        using ConnectedElement<Ieeest>::monitor_;
        using ConnectedElement<Ieeest>::signals_;

      public:
        using ScalarT    = typename ConnectedElement<Ieeest>::ScalarT;
        using IdxT       = typename ConnectedElement<Ieeest>::IdxT;
        using RealT      = typename ConnectedElement<Ieeest>::RealT;
        using ModelDataT = IeeestData<RealT, IdxT>;
        using MonitorT   = typename ConnectedElement<Ieeest>::MonitorT;

        Ieeest();
        Ieeest(const ModelDataT& data);
        ~Ieeest();

        int setGridKitComponentID(IdxT) override final;
        int allocate() override final;
        int verify() const override final;
        int initialize() override final;
        int tagDifferentiable() override final;
        int evaluateResidual() override final;
        int evaluateJacobian() override final;

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
        RealT safe_inv_T2_{1};
        RealT use_T4_block_{1};
        RealT bypass_T4_block_{0};
        RealT safe_inv_T4_{1};
        RealT use_T6_block_{1};
        RealT bypass_T6_block_{0};
        RealT safe_inv_T6_{1};
        RealT use_cutout_{1};
        RealT bypass_cutout_{0};

        void initializeParameters(const ModelDataT& data);
        void initializeMonitor();

        std::vector<ScalarT> ws_;
        std::vector<IdxT>    ws_indices_;
      };

    } // namespace Stabilizer
  } // namespace PhasorDynamics
} // namespace GridKit
