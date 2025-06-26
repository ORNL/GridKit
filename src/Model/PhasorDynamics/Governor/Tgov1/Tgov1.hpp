/**
 * @file Tgov1.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @author Adam Birchfield (abirchfield@tamu.edu)
 * @brief Declaration of a Turbine Governor Model (IEEET1).
 *
 */

#pragma once

#include <Model/PhasorDynamics/Component.hpp>

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
    class BusControl;

  } // namespace PhasorDynamics
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Governor
    {

      template <class ScalarT, typename IdxT>
      class Tgov1 : public Component<ScalarT, IdxT>
      {
        using Component<ScalarT, IdxT>::alpha_;
        using Component<ScalarT, IdxT>::f_;
        using Component<ScalarT, IdxT>::fB_;
        using Component<ScalarT, IdxT>::g_;
        using Component<ScalarT, IdxT>::gB_;
        using Component<ScalarT, IdxT>::nnz_;
        using Component<ScalarT, IdxT>::param_;
        using Component<ScalarT, IdxT>::size_;
        using Component<ScalarT, IdxT>::size_param_;
        using Component<ScalarT, IdxT>::tag_;
        using Component<ScalarT, IdxT>::time_;
        using Component<ScalarT, IdxT>::y_;
        using Component<ScalarT, IdxT>::yB_;
        using Component<ScalarT, IdxT>::yp_;
        using Component<ScalarT, IdxT>::ypB_;

        using real_type       = typename Component<ScalarT, IdxT>::real_type;
        using model_data_type = Tgov1Data<real_type, IdxT>;
        using bus_type        = BusControl<ScalarT, IdxT>;

      public:
        Tgov1(const model_data_type& data);
        Tgov1();
        ~Tgov1() = default;

        int allocate() override;
        int initialize() override;
        int tagDifferentiable() override;
        int evaluateResidual() override;

        // Still to be implemented
        int evaluateJacobian() override;
        int evaluateIntegrand() override;
        int initializeAdjoint() override;
        int evaluateAdjointResidual() override;
        int evaluateAdjointIntegrand() override;

        void updateTime(real_type /* t */, real_type /* a */) override
        {
        }

        // Setters for input signals
        void set_speed_signal(bus_type* signal);
        void set_pmech_signal(bus_type* signal);

      private:
        // Input Signals
        bus_type* speed_signal_;

        // Output Signals
        bus_type* pmech_signal_;

        // Input parameters
        real_type R_;
        real_type Pvmin_;
        real_type Pvmax_;
        real_type T1_;
        real_type T2_;
        real_type T3_;
        real_type Dt_;

        // Input States (which can be parameters)
        ScalarT pref_;

        // Scale of Sigmoid function (temporary local implementation)
        const ScalarT mu_ = 4000.0;

        // Activation function (sigmoid approximation)
        ScalarT sigmoid(ScalarT x);

        // Indicator of Valve limit states
        ScalarT indicator_low(ScalarT x, ScalarT f);
        ScalarT indicator_high(ScalarT x, ScalarT f);
        ScalarT indicator(ScalarT x, ScalarT f);
      };

    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit
