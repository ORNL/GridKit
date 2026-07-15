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
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>

// Forward declarations
namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Governor
    {
      template <typename real_type, typename index_type>
      struct Tgov1Data;
    } // namespace Governor

    template <typename scalar_type, typename index_type>
    class Genrou;

    template <typename scalar_type, typename index_type>
    class SignalNode;

  } // namespace PhasorDynamics
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Governor
    {
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

      template <typename scalar_type, typename index_type>
      class Tgov1 : public Component<scalar_type, index_type>
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
        using Component<scalar_type, index_type>::va_system_base_;
        using Component<scalar_type, index_type>::variable_indices_;
        using Component<scalar_type, index_type>::residual_indices_;
        using Component<scalar_type, index_type>::allocated_;

      public:
        using ScalarT    = scalar_type;
        using IdxT       = index_type;
        using RealT      = typename Component<ScalarT, IdxT>::RealT;
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
        int setAbsoluteTolerance(RealT) override final;
        int evaluateResidual() override final;

        // Still to be implemented
        int evaluateJacobian() override final;

        /// Get the `ComponentSignals` from this `Tgov1`
        auto getSignals()
            -> ComponentSignals<ScalarT,
                                IdxT,
                                Tgov1InternalVariables,
                                Tgov1ExternalVariables>&
        {
          return signals_;
        }

      public:
        __attribute__((always_inline)) inline int evaluateInternalResidual(
            const ScalarT*, const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);

      private:
        // Input parameters
        RealT Trate_{0};
        RealT R_{0};
        RealT Pvmin_{0};
        RealT Pvmax_{0};
        RealT T1_{0};
        RealT T2_{0};
        RealT T3_{0};
        RealT Dt_{0};

        // Derived parameters
        RealT va_component_base_{0};

        // Input States (which can be parameters)
        ScalarT pref_{0};

        /// Component signal extension
        ComponentSignals<ScalarT, IdxT, Tgov1InternalVariables, Tgov1ExternalVariables> signals_;

        // Parameter initialization function
        void    initializeParameters(const ModelDataT& data);
        void    setDerivedParams();
        ScalarT toComponentBase(ScalarT value) const;
        ScalarT toSystemBase(ScalarT value) const;

        /* Local copies of signal variables */
        std::vector<ScalarT> ws_;
        std::vector<IdxT>    ws_indices_;
      };

    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit
