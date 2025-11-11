/* Bus Fault Component - Adam Birchfield */
#pragma once

#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>

// Forward declaration of BusData structure
namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename RealT, typename IdxT>
    struct BusFaultData;
  }
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <class ScalarT, typename IdxT>
    class BusFault : public Component<ScalarT, IdxT>
    {
      using Component<ScalarT, IdxT>::gridkit_component_id_;
      using Component<ScalarT, IdxT>::alpha_;
      using Component<ScalarT, IdxT>::nnz_;
      using Component<ScalarT, IdxT>::size_;
      using Component<ScalarT, IdxT>::tag_;
      using Component<ScalarT, IdxT>::time_;
      using Component<ScalarT, IdxT>::y_;
      using Component<ScalarT, IdxT>::yp_;
      using Component<ScalarT, IdxT>::wb_;
      using Component<ScalarT, IdxT>::h_;

      using bus_type = BusBase<ScalarT, IdxT>;
      using RealT    = typename Component<ScalarT, IdxT>::RealT;
      using DataT    = BusFaultData<RealT, IdxT>;

    public:
      BusFault(bus_type* bus);
      BusFault(bus_type* bus, RealT R, RealT X, int status);
      BusFault(bus_type* bus, const DataT& data);
      ~BusFault() = default;

      int setGridKitComponentID(IdxT) override;
      int allocate() override;
      int initialize() override;
      int tagDifferentiable() override;
      int evaluateResidual() override;
      int evaluateJacobian() override;

      int verify() const override
      {
        return 0;
      }

      void updateTime(RealT /* t */, RealT /* a */) override
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
      bus_type* bus_;
      RealT     R_{0.0};
      RealT     X_{0.0};
      bool      status_{false};
      IdxT      bus_id_{0};

      /* Derivied parameters */
      RealT B_;
      RealT G_;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
