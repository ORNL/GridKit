/**
 * @file GenClassical.hpp
 * @author Abdourahman Barry (abdourahman@vt.edu)
 * @author Slaven Peles (peless@ornl.gov)
 *
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ConnectedElement.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GenClassical/GenClassicalData.hpp>

// Forward declarations.
namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename ScalarP, typename IdxP>
    class BusBase;

    template <typename RealP, typename IdxP>
    struct GenClassicalData;
  } // namespace PhasorDynamics
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename, typename>
    class GenClassical;

    enum class GenClassicalInternalVariables : size_t
    {
      MAXIMUM
    };

    enum class GenClassicalExternalVariables : size_t
    {
      MAXIMUM
    };

    template <typename ScalarP, typename IdxP>
    struct ConnectedElementTraits<GenClassical<ScalarP, IdxP>>
    {
      using GenClassicalT = GenClassical<ScalarP, IdxP>;

      using ElementT           = GenClassicalT;
      using ScalarT            = ScalarP;
      using IdxT               = IdxP;
      using RealT              = typename ScalarTraits<ScalarT>::RealT;
      using ModelDataT         = GenClassicalData<RealT, IdxT>;
      using InternalVariablesT = GenClassicalInternalVariables;
      using ExternalVariablesT = GenClassicalExternalVariables;
      using InterfaceT         = Component<ScalarT, IdxT>;
    };

    template <typename ScalarP, typename IdxP>
    class GenClassical : public ConnectedElement<GenClassical<ScalarP, IdxP>>
    {
      using ConnectedElement<GenClassical>::gridkit_component_id_;
      using ConnectedElement<GenClassical>::alpha_;
      using ConnectedElement<GenClassical>::f_;
      using ConnectedElement<GenClassical>::nnz_;
      using ConnectedElement<GenClassical>::size_;
      using ConnectedElement<GenClassical>::tag_;
      using ConnectedElement<GenClassical>::time_;
      using ConnectedElement<GenClassical>::y_;
      using ConnectedElement<GenClassical>::yp_;
      using ConnectedElement<GenClassical>::wb_;
      using ConnectedElement<GenClassical>::h_;
      using ConnectedElement<GenClassical>::J_;
      using ConnectedElement<GenClassical>::J_rows_buffer_;
      using ConnectedElement<GenClassical>::J_cols_buffer_;
      using ConnectedElement<GenClassical>::J_vals_buffer_;
      using ConnectedElement<GenClassical>::mva_system_base_;
      using ConnectedElement<GenClassical>::variable_indices_;
      using ConnectedElement<GenClassical>::residual_indices_;
      using ConnectedElement<GenClassical>::monitor_;
      using ConnectedElement<GenClassical>::signals_;

    public:
      using ScalarT    = typename ConnectedElement<GenClassical>::ScalarT;
      using IdxT       = typename ConnectedElement<GenClassical>::IdxT;
      using RealT      = typename ConnectedElement<GenClassical>::RealT;
      using BusT       = BusBase<ScalarT, IdxT>;
      using MonitorT   = typename ConnectedElement<GenClassical>::MonitorT;
      using ModelDataT = GenClassicalData<RealT, IdxT>;

      GenClassical(BusT* bus, int unit_id);
      GenClassical(BusT* bus,
                   int   unit_id,
                   RealT p0,
                   RealT q0,
                   RealT H,
                   RealT D,
                   RealT Ra,
                   RealT Xdp);
      GenClassical(BusT* bus, const ModelDataT& data);
      ~GenClassical();

      int setGridKitComponentID(IdxT) override final;
      int allocate() override final;
      int initialize() override final;
      int tagDifferentiable() override final;
      int evaluateResidual() override final;

      int verify() const override final
      {
        return 0;
      }

      // Still to be implemented
      int evaluateJacobian() override final;

      void setPmech(RealT pmech)
      {
        pmech_set_ = pmech;
      }

      void setEp(RealT ep)
      {
        ep_set_ = ep;
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
      BusT* bus_;
      IdxT  bus_id_{0};
      int   unit_id_; //< @todo this should be removed

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
    };

  } // namespace PhasorDynamics
} // namespace GridKit
