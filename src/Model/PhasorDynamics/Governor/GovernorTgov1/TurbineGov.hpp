/**
 * @file TurbineGov.hpp
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
    template <typename RealT, typename IdxT>
    struct TurbineGovData;
    template <class ScalarT, typename IdxT>
    class Genrou;
  } // namespace PhasorDynamics
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {

    template <class ScalarT, typename IdxT>
    class TurbineGov : public Component<ScalarT, IdxT>
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

      using machine_type    = Genrou<ScalarT, IdxT>;
      using real_type       = typename Component<ScalarT, IdxT>::real_type;
      using model_data_type = TurbineGovData<real_type, IdxT>;

    public:
      TurbineGov(machine_type* machine, const model_data_type& data);
      TurbineGov(machine_type* machine);
      ~TurbineGov() = default;

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

      // Read Access to Pmech
      ScalarT& Pmech();

      // Activation function (sigmoid approximation)
      ScalarT sigmoid(ScalarT x);

      // Indicator of Valve limit states
      ScalarT indicator_low(ScalarT x, ScalarT f);
      ScalarT indicator_high(ScalarT x, ScalarT f);
      ScalarT indicator(ScalarT x, ScalarT f);

    private:
    
      // Associated Machine Model
      machine_type* machine_;

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
    };

  } // namespace PhasorDynamics
} // namespace GridKit
