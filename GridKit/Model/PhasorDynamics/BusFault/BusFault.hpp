/* Bus Fault Component - Adam Birchfield */
#pragma once

#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ConnectedElement.hpp>

// Forward declaration of BusData structure
namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename RealP, typename IdxP>
    struct BusFaultData;
  }
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename, typename>
    class BusFault;

    enum class BusFaultInternalVariables : size_t
    {
      MAXIMUM
    };

    enum class BusFaultExternalVariables : size_t
    {
      MAXIMUM
    };

    template <typename ScalarP, typename IdxP>
    struct ConnectedElementTraits<BusFault<ScalarP, IdxP>>
    {
      using BusFaultT = BusFault<ScalarP, IdxP>;

      using ElementT           = BusFaultT;
      using ScalarT            = ScalarP;
      using IdxT               = IdxP;
      using RealT              = typename ScalarTraits<ScalarT>::RealT;
      using ModelDataT         = BusFaultData<RealT, IdxT>;
      using InternalVariablesT = BusFaultInternalVariables;
      using ExternalVariablesT = BusFaultExternalVariables;
      using InterfaceT         = Component<ScalarT, IdxT>;
    };

    template <typename ScalarP, typename IdxP>
    class BusFault : public ConnectedElement<BusFault<ScalarP, IdxP>>
    {
      using ConnectedElement<BusFault>::gridkit_component_id_;
      using ConnectedElement<BusFault>::alpha_;
      using ConnectedElement<BusFault>::nnz_;
      using ConnectedElement<BusFault>::size_;
      using ConnectedElement<BusFault>::tag_;
      using ConnectedElement<BusFault>::time_;
      using ConnectedElement<BusFault>::y_;
      using ConnectedElement<BusFault>::yp_;
      using ConnectedElement<BusFault>::wb_;
      using ConnectedElement<BusFault>::h_;
      using ConnectedElement<BusFault>::J_rows_buffer_;
      using ConnectedElement<BusFault>::J_cols_buffer_;
      using ConnectedElement<BusFault>::J_vals_buffer_;
      using ConnectedElement<BusFault>::monitor_;

    public:
      using ScalarT    = typename ConnectedElement<BusFault>::ScalarT;
      using IdxT       = typename ConnectedElement<BusFault>::IdxT;
      using RealT      = typename ConnectedElement<BusFault>::RealT;
      using BusT       = BusBase<ScalarT, IdxT>;
      using ModelDataT = BusFaultData<RealT, IdxT>;
      using MonitorT   = typename ConnectedElement<BusFault>::MonitorT;

      BusFault(BusT* bus);
      BusFault(BusT* bus, RealT R, RealT X, int status);
      BusFault(BusT* bus, const ModelDataT& data);
      ~BusFault();

      int setGridKitComponentID(IdxT) override final;
      int allocate() override final;
      int initialize() override final;
      int tagDifferentiable() override final;
      int evaluateResidual() override final;
      int evaluateJacobian() override final;

      int verify() const override final
      {
        return 0;
      }

      void updateTime(RealT /* t */, RealT /* a */) override final
      {
      }

    public:
      void setR(RealT R)
      {
        R_ = R;
        setDerivedParams();
      }

      void setX(RealT X)
      {
        X_ = X;
        setDerivedParams();
      }

      void setStatus(bool status)
      {
        status_ = status;
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
      __attribute__((always_inline)) inline int evaluateBusResidual(ScalarT*, ScalarT*, ScalarT*, ScalarT*);

    private:
      BusT* bus_;
      RealT R_{0.0};
      RealT X_{0.0};
      bool  status_{false};
      IdxT  bus_id_{0};

      /* Derivied parameters */
      RealT B_;
      RealT G_;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
