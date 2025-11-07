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
      template <typename RealT, typename IdxT>
      struct Ieeet1Data;
    } // namespace Exciter

    template <class ScalarT, typename IdxT>
    class BusBase;

    template <class ScalarT, typename IdxT>
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
        MAXIMUM,
      };

      template <class ScalarT, typename IdxT>
      class Ieeet1 : public Component<ScalarT, IdxT>
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

      public:
        using RealT           = typename Component<ScalarT, IdxT>::RealT;
        using model_data_type = Ieeet1Data<RealT, IdxT>;
        using signal_type     = SignalNode<ScalarT, IdxT>;
        using bus_type        = BusBase<ScalarT, IdxT>;

        Ieeet1(bus_type* bus);
        Ieeet1(signal_type*           efd_signal,
               signal_type*           speed_signal,
               bus_type*              bus,
               const model_data_type& data);
        Ieeet1(bus_type*              bus,
               const model_data_type& data);
        ~Ieeet1() = default;

        int setGridKitComponentID(IdxT) override;
        int allocate() override;
        int verify() const override;
        int initialize() override;
        int tagDifferentiable() override;
        int evaluateResidual() override;
        int evaluateJacobian() override;

        void updateTime(RealT /* t */, RealT /* a */) override
        {
        }

        /// Get the `ComponentSignals` from this `Ieeet1`
        auto getSignals()
            -> ComponentSignals<ScalarT,
                                IdxT,
                                Ieeet1InternalVariables,
                                Ieeet1ExternalVariables>&
        {
          return signals_;
        }

        bool monitoring() const override;
        void printMonitoredVariables(std::ostream& os) const override;

      private:
        // Signal pointers
        signal_type* efd_signal_;
        signal_type* speed_signal_;
        bus_type*    bus_;

        // Model Input parameters
        RealT Tr_;      ///< Time constant for voltage sensing
        RealT Ka_;      ///< Coefficient for voltage regulation
        RealT Ta_;      ///< Time constant for voltage regulation
        RealT Ke_;      ///< Coefficient for excitation system
        RealT Te_;      ///< Time constant for excitation system
        RealT Kf_;      ///< Coefficient for feedback
        RealT Tf_;      ///< Time constant for feedback
        RealT Vrmin_;   ///< LL to voltage regulation
        RealT Vrmax_;   ///< HH to voltage regulation
        RealT E1_;      ///< Saturation parameter
        RealT E2_;      ///< Saturation parameter
        RealT Se1_;     ///< Saturation parameter
        RealT Se2_;     ///< Saturation parameter
        RealT Ispdlim_; ///< Speed limit flag indicator

        // Model Derived parameters
        // TODO -> Need to be solved for in instantiation!
        RealT SA_{0};
        RealT SB_{0};

        // External Variables that don't have models yet.
        // They are constants until then.
        ScalarT vref_{0}; // (Setpoint voltage, can be different from terminal voltage)
        ScalarT vUEL_{0};
        ScalarT vOEL_{0};
        ScalarT vS_{0};
        ScalarT Ec_{0}; // "Compensated" terminal measurment

        /// Component signal extension
        ComponentSignals<ScalarT, IdxT, Ieeet1InternalVariables, Ieeet1ExternalVariables> signals_;

        // Scale of Sigmoid function
        // (temporary local implementation)
        // This value gave higher precision.
        const RealT mu_{400000.0};

        Model::VariableMonitor<Ieeet1, Ieeet1Data> monitor_;

        // Activation function (sigmoid approximation)
        ScalarT sigmoid(ScalarT x);
        ScalarT indicator(ScalarT x, ScalarT f);

        // Parameter initialization function
        void initModelParams(const model_data_type& data);

        void initializeMonitor();
      };

    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
