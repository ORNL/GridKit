/**
 * @file GenClassical.hpp
 * @author Abdourahman Barry (abdourahman@vt.edu)
 * @author Slaven Peles (peless@ornl.gov)
 *
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GenClassical/GenClassicalData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

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

      using bus_type = BusBase<ScalarT, IdxT>;
      using RealT    = typename Component<ScalarT, IdxT>::RealT;
      using DataT    = GenClassicalData<RealT, IdxT>;

    public:
      GenClassical(bus_type* bus, int unit_id);
      GenClassical(bus_type* bus,
                   int       unit_id,
                   RealT     p0,
                   RealT     q0,
                   RealT     H,
                   RealT     D,
                   RealT     Ra,
                   RealT     Xdp);
      GenClassical(bus_type* bus, const DataT& data);
      ~GenClassical() = default;

      int setGridKitComponentID(IdxT) override;
      int allocate() override;
      int initialize() override;
      int tagDifferentiable() override;
      int evaluateResidual() override;

      int verify() const override
      {
        return 0;
      }

      // Still to be implemented
      int evaluateJacobian() override;

      void updateTime(RealT /* t */, RealT /* a */) override
      {
      }

      void setPmech(RealT pmech)
      {
        pmech_set_ = pmech;
      }

      void setEp(RealT ep)
      {
        ep_set_ = ep;
      }

      const Model::VariableMonitorBase* getMonitor() const override
      {
        return &monitor_;
      }

    private:
      void initializeMonitor();
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
      RealT p0_{0.0};
      RealT q0_{0.0};

      /* Input parameters */
      RealT H_{0.0};
      RealT D_{0.0};
      RealT Ra_{0.0};
      RealT Xdp_{0.0};
      RealT mva_base_{100.0};

      /* Derivied parameters */
      RealT G_;
      RealT B_;

      /* Setpoints for control variables (determined at initialization) */
      ScalarT pmech_set_;
      ScalarT ep_set_;

      /// Variable monitor
      Model::VariableMonitor<GenClassical, GenClassicalData> monitor_;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
