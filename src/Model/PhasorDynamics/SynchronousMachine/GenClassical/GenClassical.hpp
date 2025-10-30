/**
 * @file GenClassical.hpp
 * @author Abdourahman Barry (abdourahman@vt.edu)
 * @author Slaven Peles (peless@ornl.gov)
 *
 */

#pragma once

#include <Model/PhasorDynamics/Component.hpp>
#include <Model/PhasorDynamics/ComponentSignals.hpp>

// Forward declarations.
namespace GridKit
{
  namespace PhasorDynamics
  {
    template <class ScalarT, typename IdxT>
    class BusBase;

    template <typename RealT, typename IdxT>
    struct GenClassicalData;
  } // namespace PhasorDynamics
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {

    template <class ScalarT, typename IdxT>
    class GenClassical : public Component<ScalarT, IdxT>
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
      using Component<ScalarT, IdxT>::wb_;
      using Component<ScalarT, IdxT>::h_;
      using Component<ScalarT, IdxT>::J_;
      using Component<ScalarT, IdxT>::mva_system_base_;

      using bus_type  = BusBase<ScalarT, IdxT>;
      using real_type = typename Component<ScalarT, IdxT>::real_type;
      using DataT     = GenClassicalData<real_type, IdxT>;

    public:
      GenClassical(bus_type* bus, int unit_id);
      GenClassical(bus_type* bus,
                   int       unit_id,
                   real_type p0,
                   real_type q0,
                   real_type H,
                   real_type D,
                   real_type Ra,
                   real_type Xdp);
      GenClassical(bus_type* bus, const DataT& data);
      ~GenClassical() = default;

      int setGridKitComponentID(IdxT) override;
      int allocate() override;
      int initialize() override;
      int tagDifferentiable() override;
      int evaluateResidual() override;

      // Still to be implemented
      int evaluateJacobian() override;

      void updateTime(real_type /* t */, real_type /* a */) override
      {
      }

      void setPmech(real_type pmech)
      {
        pmech_set_ = pmech;
      }

      void setEp(real_type ep)
      {
        ep_set_ = ep;
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

    public:
      __attribute__((always_inline)) inline int evaluateInternalResidual(ScalarT*, ScalarT*, ScalarT*, ScalarT*);
      __attribute__((always_inline)) inline int evaluateBusResidual(ScalarT*, ScalarT*, ScalarT*, ScalarT*);

    private:
      /* Identification */
      bus_type* bus_;
      IdxT      bus_id_{0};
      int       unit_id_; //< @todo this should be removed

      /* Initial terminal conditions */
      real_type p0_{0.0};
      real_type q0_{0.0};

      /* Input parameters */
      real_type H_{0.0};
      real_type D_{0.0};
      real_type Ra_{0.0};
      real_type Xdp_{0.0};
      real_type mva_base_{100.0};

      /* Derivied parameters */
      real_type G_;
      real_type B_;

      /* Setpoints for control variables (determined at initialization) */
      ScalarT pmech_set_;
      ScalarT ep_set_;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
