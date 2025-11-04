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

      using bus_type  = BusBase<ScalarT, IdxT>;
      using real_type = typename Component<ScalarT, IdxT>::real_type;
      using DataT     = BusFaultData<real_type, IdxT>;

    public:
      BusFault(bus_type* bus);
      BusFault(bus_type* bus, real_type R, real_type X, int status);
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

      void updateTime(real_type /* t */, real_type /* a */) override
      {
      }

    public:
      void setR(real_type R)
      {
        R_ = R;
        setDerivedParams();
      }

      void setX(real_type X)
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
      real_type R_{0.0};
      real_type X_{0.0};
      bool      status_{false};
      IdxT      bus_id_{0};

      /* Derivied parameters */
      real_type B_;
      real_type G_;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
