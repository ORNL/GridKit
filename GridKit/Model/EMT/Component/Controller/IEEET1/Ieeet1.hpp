/**
 * @file   Ieeet1.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @author Adam Birchfield (abirchfield@tamu.edu)
 *
 * @brief Declaration of a IEEET1 Controller Model.
 *
 */

#pragma once

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Controller/IEEET1/Ieeet1Data.hpp>
#include <GridKit/Model/EMT/ComponentSignals.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

// Forward declarations
namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      template <typename real_type, typename index_type>
      struct Ieeet1Data;
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
      /// Internal variables of a `Ieeet1`
      enum class Ieeet1InternalVariables : size_t
      {
        VTS,  ///< Sensed term voltage
        VR,   ///< Voltage regulation
        EFDP, ///< Efd (pre multiplication)
        VFX,  ///< Exciter feedback
        VTR,  ///< Terminal voltage error
        VF,   ///< Feedback voltage
        VE,   ///< Excitation control voltage
        EFD,  ///< Efd
        KSAT, ///< \f$E_{\mathrm{fd}}'S(E_{\mathrm{fd}}')\f$ Scaled-quadratic saturation contribution
        MAXIMUM,
      };

      /// External variables of a `Ieeet1`
      enum class Ieeet1ExternalVariables : size_t
      {
        VA,    ///< Phase-a terminal voltage in volts
        VB,    ///< Phase-b terminal voltage in volts
        VC,    ///< Phase-c terminal voltage in volts
        OMEGA, ///< Machine rotor speed in per unit
        VREF,  ///< Voltage reference
        VS,    ///< Stabilizer output signal
        VUEL,  ///< Under-excitation limiter signal
        VOEL,  ///< Over-excitation limiter signal
        MAXIMUM,
      };

      template <typename scalar_type, typename index_type>
      class Ieeet1 : public Component<scalar_type, index_type>
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
        using ModelDataT = Ieeet1Data<RealT, IdxT>;
        using SignalT    = Signal<ScalarT, IdxT>;
        using MonitorT   = Model::VariableMonitor<Ieeet1, Ieeet1Data>;

        Ieeet1();
        explicit Ieeet1(const ModelDataT& data);
        ~Ieeet1();

        int setGridKitComponentID(IdxT) override final;
        int allocate() override final;
        int verify() const override final;

        int initializationOrder() const noexcept override final
        {
          return 2;
        }

        int initialize() override final;
        int tagDifferentiable() override final;
        int setAbsoluteTolerance(RealT rel_tol) override final;
        int evaluateInternalResidual() override final;
        int evaluateResidual() override final;
        int evaluateJacobian() override final;

        /// Get the `ComponentSignals` from this `Ieeet1`
        auto getSignals()
            -> ComponentSignals<ScalarT,
                                IdxT,
                                Ieeet1InternalVariables,
                                Ieeet1ExternalVariables>&
        {
          return signals_;
        }

        const Model::VariableMonitorBase* getMonitor() const override;

        __attribute__((always_inline)) inline int evaluateInternalResidual(
            const ScalarT*, const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);

      private:
        static constexpr RealT TIME_CONSTANT_MINIMUM = static_cast<RealT>(1.0e-3);
        static void            logTimeConstantWarning();

        // Model input parameters
        RealT V_{0.0};       ///< Rated line-to-line RMS voltage in volts
        RealT Tr_{0.0};      ///< Time constant for voltage sensing
        RealT Ka_{50.0};     ///< Coefficient for voltage regulation
        RealT Ta_{0.04};     ///< Time constant for voltage regulation
        RealT Ke_{-0.06};    ///< Coefficient for excitation system
        RealT Te_{0.6};      ///< Time constant for excitation system
        RealT Kf_{0.09};     ///< Coefficient for feedback
        RealT Tf_{1.46};     ///< Time constant for feedback
        RealT Vrmin_{-1.0};  ///< LL to voltage regulation
        RealT Vrmax_{1.0};   ///< HH to voltage regulation
        RealT E1_{2.8};      ///< Saturation parameter
        RealT E2_{3.73};     ///< Saturation parameter
        RealT Se1_{0.04};    ///< Saturation parameter
        RealT Se2_{0.33};    ///< Saturation parameter
        RealT Ispdlim_{0.0}; ///< Speed limit flag indicator

        // Model Derived parameters
        RealT Ke_eff_{Ke_};

        // Saturation coefficients derived from E1, E2, Se1, and Se2.
        RealT SA_{0};
        RealT SB_{0};

        // Runtime limiter masks keep the summing junction Enzyme sparse-solvable.
        RealT uel_on_{0};
        RealT oel_on_{0};

        ScalarT omega_set_{1};
        ScalarT vref_set_{0};
        ScalarT vs_set_{0};
        ScalarT vuel_set_{0};
        ScalarT voel_set_{0};

        /// Component signal extension
        ComponentSignals<ScalarT, IdxT, Ieeet1InternalVariables, Ieeet1ExternalVariables> signals_;

        /// Variable monitor
        std::unique_ptr<MonitorT> monitor_;

        // Parameter initialization function
        void initModelParams(const ModelDataT& data);
        void setDerivedParameters();

        /// Associate variable getter functions with enum values
        void initializeMonitor();

        __attribute__((always_inline)) inline ScalarT voltageMagnitude(const ScalarT* voltage) const
        {
          const ScalarT norm_squared = voltage[0] * voltage[0] + voltage[1] * voltage[1]
                                       + voltage[2] * voltage[2];
          if (norm_squared == ZERO<RealT>)
            return ScalarT{ZERO<RealT>};
          return std::sqrt(norm_squared) / V_;
        }

      protected:
        void gatherExternalVariables() override;
      };

    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
