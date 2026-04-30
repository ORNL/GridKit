/**
 * @file Tgov1.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @author Adam Birchfield (abirchfield@tamu.edu)
 * @author Wiktoria Zielinska (zielinskawa@ORNL.gov)
 * @brief Declaration of a Turbine Governor Model (IEEET1).
 *
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ConnectedElement.hpp>

// Forward declarations
namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Governor
    {
      template <typename RealP, typename IdxP>
      struct Tgov1Data;
    } // namespace Governor

    template <typename ScalarP, typename IdxP>
    class Genrou;

    template <typename ScalarP, typename IdxP>
    class SignalNode;

  } // namespace PhasorDynamics
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Governor
    {
      template <typename, typename>
      class Tgov1;

      /// Internal variables of a `Tgov1`
      enum class Tgov1InternalVariables : size_t
      {
        PTX, ///< $P_{tx}$
        PV,  ///< $P_v$
        PM,  ///< $P_m$
        MAXIMUM,
      };

      /// External variables of a `Tgov1`
      enum class Tgov1ExternalVariables : size_t
      {
        DELTAOMEGA, ///< $\Delta_\omega$
        PREF,       ///< $P_{ref}$
        MAXIMUM,
      };
    } // namespace Governor

    template <typename ScalarP, typename IdxP>
    struct ConnectedElementTraits<Governor::Tgov1<ScalarP, IdxP>>
    {
      using Tgov1T = Governor::Tgov1<ScalarP, IdxP>;

      using ElementT           = Tgov1T;
      using ScalarT            = ScalarP;
      using IdxT               = IdxP;
      using RealT              = typename ScalarTraits<ScalarT>::RealT;
      using ModelDataT         = Governor::Tgov1Data<RealT, IdxT>;
      using InternalVariablesT = Governor::Tgov1InternalVariables;
      using ExternalVariablesT = Governor::Tgov1ExternalVariables;
      using InterfaceT         = Component<ScalarT, IdxT>;
    };

    namespace Governor
    {
      template <typename ScalarP, typename IdxP>
      class Tgov1 : public ConnectedElement<Tgov1<ScalarP, IdxP>>
      {
        using ConnectedElement<Tgov1>::gridkit_component_id_;
        using ConnectedElement<Tgov1>::alpha_;
        using ConnectedElement<Tgov1>::f_;
        using ConnectedElement<Tgov1>::nnz_;
        using ConnectedElement<Tgov1>::size_;
        using ConnectedElement<Tgov1>::tag_;
        using ConnectedElement<Tgov1>::time_;
        using ConnectedElement<Tgov1>::y_;
        using ConnectedElement<Tgov1>::yp_;
        using ConnectedElement<Tgov1>::wb_;
        using ConnectedElement<Tgov1>::h_;
        using ConnectedElement<Tgov1>::J_;
        using ConnectedElement<Tgov1>::J_rows_buffer_;
        using ConnectedElement<Tgov1>::J_cols_buffer_;
        using ConnectedElement<Tgov1>::J_vals_buffer_;
        using ConnectedElement<Tgov1>::variable_indices_;
        using ConnectedElement<Tgov1>::residual_indices_;
        using ConnectedElement<Tgov1>::signals_;

      public:
        using ScalarT    = typename ConnectedElement<Tgov1>::ScalarT;
        using IdxT       = typename ConnectedElement<Tgov1>::IdxT;
        using RealT      = typename ConnectedElement<Tgov1>::RealT;
        using ModelDataT = Tgov1Data<RealT, IdxT>;
        using SignalT    = SignalNode<ScalarT, IdxT>;

        Tgov1();
        Tgov1(SignalT*, SignalT*);
        Tgov1(const ModelDataT&);
        ~Tgov1() = default;

        int setGridKitComponentID(IdxT) override final;
        int allocate() override final;
        int verify() const override final;
        int initialize() override final;
        int tagDifferentiable() override final;
        int evaluateResidual() override final;

        // Still to be implemented
        int evaluateJacobian() override final;

      public:
        __attribute__((always_inline)) inline int evaluateInternalResidual(ScalarT*, ScalarT*, ScalarT*, ScalarT*, ScalarT*);

      private:
        // Input parameters
        RealT R_{0};
        RealT Pvmin_{0};
        RealT Pvmax_{0};
        RealT T1_{0};
        RealT T2_{0};
        RealT T3_{0};
        RealT Dt_{0};

        // Input States (which can be parameters)
        ScalarT pref_{0};

        // Parameter initialization function
        void initializeParameters(const ModelDataT& data);

        /* Local copies of signal variables */
        std::vector<ScalarT> ws_;
        std::vector<IdxT>    ws_indices_;
      };

    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit
