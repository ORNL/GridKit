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
#include <GridKit/Model/PhasorDynamics/ConnectedElement.hpp>

// Forward declarations
namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      template <typename RealP, typename IdxP>
      struct Ieeet1Data;
    } // namespace Exciter

    template <typename ScalarP, typename IdxP>
    class BusBase;

    template <typename ScalarP, typename IdxP>
    class SignalNode;

  } // namespace PhasorDynamics
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      template <typename, typename>
      class Ieeet1;

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
    } // namespace Exciter

    template <typename ScalarP, typename IdxP>
    struct ConnectedElementTraits<Exciter::Ieeet1<ScalarP, IdxP>>
    {
      using Ieeet1T = Exciter::Ieeet1<ScalarP, IdxP>;

      using ElementT           = Ieeet1T;
      using ScalarT            = ScalarP;
      using IdxT               = IdxP;
      using RealT              = typename ScalarTraits<ScalarT>::RealT;
      using ModelDataT         = Exciter::Ieeet1Data<RealT, IdxT>;
      using InternalVariablesT = Exciter::Ieeet1InternalVariables;
      using ExternalVariablesT = Exciter::Ieeet1ExternalVariables;
      using InterfaceT         = Component<ScalarT, IdxT>;
    };

    namespace Exciter
    {
      template <typename ScalarP, typename IdxP>
      class Ieeet1 : public ConnectedElement<Ieeet1<ScalarP, IdxP>>
      {
        using ConnectedElement<Ieeet1>::gridkit_component_id_;
        using ConnectedElement<Ieeet1>::alpha_;
        using ConnectedElement<Ieeet1>::f_;
        using ConnectedElement<Ieeet1>::nnz_;
        using ConnectedElement<Ieeet1>::size_;
        using ConnectedElement<Ieeet1>::tag_;
        using ConnectedElement<Ieeet1>::time_;
        using ConnectedElement<Ieeet1>::y_;
        using ConnectedElement<Ieeet1>::yp_;
        using ConnectedElement<Ieeet1>::wb_;
        using ConnectedElement<Ieeet1>::J_;
        using ConnectedElement<Ieeet1>::J_rows_buffer_;
        using ConnectedElement<Ieeet1>::J_cols_buffer_;
        using ConnectedElement<Ieeet1>::J_vals_buffer_;
        using ConnectedElement<Ieeet1>::variable_indices_;
        using ConnectedElement<Ieeet1>::residual_indices_;
        using ConnectedElement<Ieeet1>::signals_;
        using ConnectedElement<Ieeet1>::monitor_;

      public:
        using ScalarT    = typename ConnectedElement<Ieeet1>::ScalarT;
        using IdxT       = typename ConnectedElement<Ieeet1>::IdxT;
        using RealT      = typename ConnectedElement<Ieeet1>::RealT;
        using ModelDataT = Ieeet1Data<RealT, IdxT>;
        using SignalT    = SignalNode<ScalarT, IdxT>;
        using BusT       = BusBase<ScalarT, IdxT>;
        using MonitorT   = typename ConnectedElement<Ieeet1>::MonitorT;

        Ieeet1(BusT* bus);

        Ieeet1(SignalT*          efd_signal,
               SignalT*          speed_signal,
               BusT*             bus,
               const ModelDataT& data);

        Ieeet1(BusT* bus, const ModelDataT& data);

        ~Ieeet1();

        int setGridKitComponentID(IdxT) override final;
        int allocate() override final;
        int verify() const override final;
        int initialize() override final;
        int tagDifferentiable() override final;
        int evaluateResidual() override final;
        int evaluateJacobian() override final;

        __attribute__((always_inline)) inline int evaluateInternalResidual(ScalarT*, ScalarT*, ScalarT*, ScalarT*, ScalarT*);

      private:
        // Signal pointers
        SignalT* efd_signal_;
        SignalT* speed_signal_;
        BusT*    bus_;

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
        ScalarT Ec_{0}; // "Compensated" terminal measurment, currently unused

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
