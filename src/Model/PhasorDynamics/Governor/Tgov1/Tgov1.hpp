/**
 * @file Tgov1.hpp
 * @author Wiktoria Zielinska (zielinskawa@ORNL.gov)
 * @author Luke Lowery (lukel@tamu.edu)
 * @author Adam Birchfield (abirchfield@tamu.edu)
 * @brief Declaration of a Turbine Governor Model (IEEET1).
 *
 */

#pragma once

#include <Model/PhasorDynamics/Component.hpp>
#include <Model/PhasorDynamics/GovernorBase.hpp>

// Forward declarations
namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Governor
    {
      template <typename RealT, typename IdxT>
      struct Tgov1Data;
    } // namespace Governor

    template <class ScalarT, typename IdxT>
    class Genrou;

    template <class ScalarT, typename IdxT>
    class SignalNode;

  } // namespace PhasorDynamics
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Governor
    {

      template <class ScalarT, typename IdxT>
      class Tgov1 : public Component<ScalarT, IdxT>, public GovernorBase<ScalarT, IdxT>
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
        using model_data_type = Tgov1Data<real_type, IdxT>;
        using signal_type     = SignalNode<ScalarT, IdxT>;
        using machine_type    = Genrou<ScalarT, IdxT>;

      public:
        Tgov1(signal_type* pmech, signal_type* omega, const model_data_type& data);
        Tgov1(signal_type* pmech, signal_type* omega);
        Tgov1(machine_type* machine, const model_data_type& data);
        ~Tgov1() = default;

        int allocate() override;
        int initialize() override;
        int tagDifferentiable() override;
        int evaluateResidual() override;

        // Still to be implemented
        int evaluateJacobian() override;

        void updateTime(real_type /* t */, real_type /* a */) override
        {
        }

        // Read Access to Pmech
        ScalarT& Pmech() override;

      private:
        // Associated Machine Model
        machine_type* machine_{nullptr};
        signal_type*  pmech_{nullptr};
        signal_type*  omega_{nullptr};

        // Input parameters
        real_type R_{0.05};
        real_type Pvmin_{0.0};
        real_type Pvmax_{1.0};
        real_type T1_{0.5};
        real_type T2_{2.5};
        real_type T3_{7.5};
        real_type Dt_{0.0};

        // Input States (which can be parameters)
        ScalarT pref_{0};

        // Scale of Sigmoid function (temporary local implementation)
        const ScalarT mu_{4000.0};

        // Activation function (sigmoid approximation)
        ScalarT sigmoid(ScalarT x);

        // Indicator of Valve limit states
        ScalarT indicator_low(ScalarT x, ScalarT f);
        ScalarT indicator_high(ScalarT x, ScalarT f);
        ScalarT indicator(ScalarT x, ScalarT f);

        void initializeParameters(const model_data_type& data);
      };

    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit
