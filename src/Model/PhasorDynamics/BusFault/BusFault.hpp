/* Bus Fault Component - Adam Birchfield */
#pragma once

#include <Model/PhasorDynamics/BusBase.hpp>
#include <Model/PhasorDynamics/Component.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <class ScalarT, typename IdxT>
    class BusFault : public Component<ScalarT, IdxT>
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

      using bus_type  = BusBase<ScalarT, IdxT>;
      using real_type = typename Component<ScalarT, IdxT>::real_type;

    public:
      BusFault(bus_type* bus);
      BusFault(bus_type* bus, real_type R, real_type X, int status);
      ~BusFault() = default;

      int allocate() override;
      int initialize() override;
      int tagDifferentiable() override;
      int evaluateResidual() override;
      int evaluateJacobian() override;
      int evaluateIntegrand() override;
      int initializeAdjoint() override;
      int evaluateAdjointResidual() override;
      int evaluateAdjointIntegrand() override;

      void updateTime(real_type /* t */, real_type /* a */) override
      {
      }

    public:
      void setR(real_type R)
      {
        R_ = R;
      }

      void setX(real_type X)
      {
        X_ = X;
      }

      void setStatus(int status)
      {
        status_ = status;
      }

    private:
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
      bus_type* bus_;
      real_type R_;
      real_type X_;
      int       status_;
      const int busID_;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
