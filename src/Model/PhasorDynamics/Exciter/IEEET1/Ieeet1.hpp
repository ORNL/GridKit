/**
 * @file   Ieeet1.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @author Adam Birchfield (abirchfield@tamu.edu)
 *
 * @brief Declaration of a IEEET1 Exciter Model.
 *
 */

#pragma once

#include <Model/PhasorDynamics/Component.hpp>

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

      template <class ScalarT, typename IdxT>
      class Ieeet1 : public Component<ScalarT, IdxT>
      {
        using Component<ScalarT, IdxT>::alpha_;
        using Component<ScalarT, IdxT>::f_;
        using Component<ScalarT, IdxT>::nnz_;
        using Component<ScalarT, IdxT>::size_;
        using Component<ScalarT, IdxT>::tag_;
        using Component<ScalarT, IdxT>::time_;
        using Component<ScalarT, IdxT>::y_;
        using Component<ScalarT, IdxT>::yp_;

        using real_type       = typename Component<ScalarT, IdxT>::real_type;
        using model_data_type = Ieeet1Data<real_type, IdxT>;
        using signal_type     = SignalNode<ScalarT, IdxT>;
        using bus_type        = BusBase<ScalarT, IdxT>;

      public:
        Ieeet1(
            signal_type*           efd_signal,
            signal_type*           speed_signal,
            bus_type*              bus,
            const model_data_type& data);
        ~Ieeet1() = default;

        int allocate() override;
        int initialize() override;
        int tagDifferentiable() override;
        int evaluateResidual() override;
        int evaluateJacobian() override;

        void updateTime(real_type /* t */, real_type /* a */) override
        {
        }

      private:
        // Signal pointers
        signal_type* efd_signal_;
        signal_type* speed_signal_;
        bus_type*    bus_;

        // Model Input parameters
        real_type Tr_;      ///< Time constant for voltage sensing
        real_type Ka_;      ///< Coefficient for voltage regulation
        real_type Ta_;      ///< Time constant for voltage regulation
        real_type Ke_;      ///< Coefficient for excitation system
        real_type Te_;      ///< Time constant for excitation system
        real_type Kf_;      ///< Coefficient for feedback
        real_type Tf_;      ///< Time constant for feedback
        real_type Vrmin_;   ///< LL to voltage regulation
        real_type Vrmax_;   ///< HH to voltage regulation
        real_type E1_;      ///< Saturation parameter
        real_type E2_;      ///< Saturation parameter
        real_type Se1_;     ///< Saturation parameter
        real_type Se2_;     ///< Saturation parameter
        real_type Ispdlim_; ///< Speed limit flag indicator

        // Model Derived parameters
        // TODO -> Need to be solved for in instantiation!
        real_type SA_{0};
        real_type SB_{0};

        // External Variables that don't have models yet.
        // They are constants until then.
        ScalarT vref_{0}; // (Setpoint voltage, can be different from terminal voltage)
        ScalarT vUEL_{0};
        ScalarT vOEL_{0};
        ScalarT vS_{0};
        ScalarT Ec_{0}; // "Compensated" terminal measurment

        // Scale of Sigmoid function
        // (temporary local implementation)
        // This value gave higher precision.
        const real_type mu_{400000.0};

        // Activation function (sigmoid approximation)
        ScalarT sigmoid(ScalarT x);
        ScalarT indicator(ScalarT x, ScalarT f);

        // Parameter initialization function
        void initModelParams(const model_data_type& data);
      };

    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
