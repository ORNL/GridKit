/**
 * @file Tgov1.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the EMT TGOV1 turbine-governor model.
 *
 */

#pragma once

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Controller/TGOV1/Tgov1Data.hpp>
#include <GridKit/Model/EMT/ComponentSignals.hpp>

// Forward declarations
namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      template <typename real_type, typename index_type>
      struct Tgov1Data;
    } // namespace Controller

    template <typename scalar_type, typename index_type>
    class Signal;

  } // namespace EMT
} // namespace GridKit

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      /// Internal variables of a `Tgov1`
      enum class Tgov1InternalVariables : size_t
      {
        PTX, ///< \f$P_t\f$ Turbine-block output
        PV,  ///< \f$P_v\f$ Valve position
        PM,  ///< \f$P_m\f$ Mechanical-power output
        MAXIMUM,
      };

      /// External variables of a `Tgov1`
      enum class Tgov1ExternalVariables : size_t
      {
        OMEGA, ///< \f$\omega_r\f$ Machine rotor speed
        PREF,  ///< \f$P_\mathrm{ref}\f$ Governor reference
        MAXIMUM,
      };

      /*!
       * @brief Implementation of the TGOV1 turbine governor.
       *
       * The model is the PhasorDynamics TGOV1 governor carried into the EMT
       * subsystem. The speed input is the machine rotor speed; the speed
       * deviation from synchronous drives the droop. Mechanical power is
       * exchanged in the machine per-unit base.
       */
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
        using Component<scalar_type, index_type>::y_ext_;
        using Component<scalar_type, index_type>::yp_ext_;
        using Component<scalar_type, index_type>::variable_indices_ext_;
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
        using ModelDataT = Tgov1Data<RealT, IdxT>;
        using SignalT    = Signal<ScalarT, IdxT>;

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
        int evaluateInternalResidual() override final;
        int evaluateResidual() override final;

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
        RealT Trate_{static_cast<RealT>(100.0)};
        RealT R_{static_cast<RealT>(0.05)};
        RealT Pvmin_{ZERO<RealT>};
        RealT Pvmax_{ONE<RealT>};
        RealT T1_{static_cast<RealT>(0.5)};
        RealT T2_{static_cast<RealT>(2.5)};
        RealT T3_{static_cast<RealT>(7.5)};
        RealT Dt_{ZERO<RealT>};

        // Input States (which can be parameters)
        ScalarT pref_set_{0};

        /// Component signal extension
        ComponentSignals<ScalarT, IdxT, Tgov1InternalVariables, Tgov1ExternalVariables> signals_;

        // Parameter initialization function
        void initializeParameters(const ModelDataT& data);
        void setDerivedParams();

        static constexpr RealT TIME_CONSTANT_MINIMUM = static_cast<RealT>(1.0e-3);
        static void            logTimeConstantWarning();

      protected:
        void gatherExternalVariables() override;
      };

    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
