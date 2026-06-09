/**
 * @file   Ieeet1.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @author Adam Birchfield (abirchfield@tamu.edu)
 *
 * @brief Declaration of a IEEET1 Exciter Model.
 *
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

// Forward declarations
namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      template <typename real_type, typename index_type>
      struct Ieeet1Data;
    } // namespace Exciter

    template <typename scalar_type, typename index_type>
    class BusBase;

    template <typename scalar_type, typename index_type>
    class SignalNode;

  } // namespace PhasorDynamics
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
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
        VE,   ///< Exciter control voltage
        EFD,  ///< Efd
        KSAT, ///< Saturation
        MAXIMUM,
      };

      /// External variables of a `Ieeet1`
      enum class Ieeet1ExternalVariables : size_t
      {
        OMEGA, ///< Generator speed deviation
        VREAL, ///< Real bus voltage
        VIMAG, ///< Imaginary bus voltage
        VS,    ///< Stabilizer output signal
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
        using Component<scalar_type, index_type>::wb_;
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
        using BusT       = BusBase<ScalarT, IdxT>;
        using ModelDataT = Ieeet1Data<RealT, IdxT>;
        using SignalT    = SignalNode<ScalarT, IdxT>;
        using MonitorT   = Model::VariableMonitor<Ieeet1, Ieeet1Data>;

        Ieeet1(BusT* bus);
        Ieeet1(BusT*             bus,
               const ModelDataT& data);
        ~Ieeet1();

        int setGridKitComponentID(IdxT) override final;
        int allocate() override final;
        int verify() const override final;
        int initialize() override final;
        int tagDifferentiable() override final;
        int setAbsoluteTolerance(RealT rel_tol) override final;
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

        FORCE_INLINE int evaluateInternalResidual(ScalarT*, ScalarT*, ScalarT*, ScalarT*, ScalarT*);

      private:
        static constexpr RealT TIME_CONSTANT_MINIMUM = static_cast<RealT>(1.0e-3);

        // Signal pointers
        BusT* bus_;

        // Model Input parameters
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
        // Saturation coefficients derived from E1, E2, Se1, and Se2.
        RealT SA_{0};
        RealT SB_{0};

        // External Variables that don't have models yet.
        // They are constants until then.
        ScalarT vref_{0}; // (Setpoint voltage, can be different from terminal voltage)
        ScalarT vUEL_{0};
        ScalarT vOEL_{0};

        /// Component signal extension
        ComponentSignals<ScalarT, IdxT, Ieeet1InternalVariables, Ieeet1ExternalVariables> signals_;

        /// Variable monitor
        std::unique_ptr<MonitorT> monitor_;

        // Parameter initialization function
        void initModelParams(const ModelDataT& data);

        /// Associate variable getter functions with enum values
        void initializeMonitor();

        /* Local copies of signal variables */
        std::vector<ScalarT> ws_;
        std::vector<IdxT>    ws_indices_;
      };

    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
