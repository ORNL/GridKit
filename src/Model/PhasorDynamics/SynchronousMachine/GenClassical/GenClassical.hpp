/**
 * @file GenClassical.cpp
 * @author Abdourahman Barry (abdourahman@vt.edu)
 *
 */

#pragma once

#include <Model/PhasorDynamics/Component.hpp>

// Forward declarations.
namespace GridKit
{
  namespace PhasorDynamics
  {
    template <class ScalarT, typename IdxT>
    class BusElectric;
  }
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {

    template <class ScalarT, typename IdxT>
    class GenClassical : public Component<ScalarT, IdxT>
    {
      using Component<ScalarT, IdxT>::alpha_;
      using Component<ScalarT, IdxT>::f_;
      using Component<ScalarT, IdxT>::fB_;
      using Component<ScalarT, IdxT>::g_;
      using Component<ScalarT, IdxT>::gB_;
      using Component<ScalarT, IdxT>::nnz_;
      using Component<ScalarT, IdxT>::param_;
      using Component<ScalarT, IdxT>::size_;
      using Component<ScalarT, IdxT>::tag_;
      using Component<ScalarT, IdxT>::time_;
      using Component<ScalarT, IdxT>::y_;
      using Component<ScalarT, IdxT>::yB_;
      using Component<ScalarT, IdxT>::yp_;
      using Component<ScalarT, IdxT>::ypB_;

      using bus_type  = BusElectric<ScalarT, IdxT>;
      using real_type = typename Component<ScalarT, IdxT>::real_type;

    public:
      GenClassical(bus_type* bus, int unit_id);
      GenClassical(bus_type* bus,
                   int       unit_id,
                   ScalarT   p0,
                   ScalarT   q0,
                   real_type H,
                   real_type D,
                   real_type Ra,
                   real_type Xdp);
      ~GenClassical() = default;

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

    private:
      void setDerivedParams();

      ScalarT& Vr()
      {
        return bus_->Vr();
      }

      ScalarT& Vi()
      {
        return bus_->Vi();
      }

      ScalarT& Ir()
      {
        return bus_->Ir();
      }

      ScalarT& Ii()
      {
        return bus_->Ii();
      }

    private:
      /* Identification */
      bus_type* bus_;
      const int busID_;
      int       unit_id_;

      /* Initial terminal conditions */
      ScalarT p0_;
      ScalarT q0_;

      /* Input parameters */
      real_type H_;
      real_type D_;
      real_type Ra_;
      real_type Xdp_;

      /* Derivied parameters */
      real_type G_;
      real_type B_;

      /* Setpoints for control variables (determined at initialization) */
      real_type pmech_set_;
      real_type ep_set_;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
