
#ifndef _BUS_PV_HPP_
#define _BUS_PV_HPP_

#include <cassert>

#include "BaseBus.hpp"
#include <Model/PowerFlow/PowerSystemData.hpp>

namespace GridKit
{
  /*!
   * @brief Implementation of a PV bus.
   *
   * Voltage _V_ and phase _theta_ are variables in PV bus model.
   * Active and reactive power, _P_ and _Q_, are residual components.
   *
   *
   */
  template <class ScalarT, typename IdxT>
  class BusPV : public BaseBus<ScalarT, IdxT>
  {
    using BaseBus<ScalarT, IdxT>::size_;
    using BaseBus<ScalarT, IdxT>::y_;
    using BaseBus<ScalarT, IdxT>::yp_;
    using BaseBus<ScalarT, IdxT>::yB_;
    using BaseBus<ScalarT, IdxT>::ypB_;
    using BaseBus<ScalarT, IdxT>::f_;
    using BaseBus<ScalarT, IdxT>::fB_;
    using BaseBus<ScalarT, IdxT>::tag_;

  public:
    using real_type = typename ModelEvaluatorImpl<ScalarT, IdxT>::real_type;
    using BusData   = GridKit::PowerSystemData::BusData<real_type, IdxT>;

    BusPV();
    BusPV(ScalarT V, ScalarT theta0);
    BusPV(BusData& data);
    virtual ~BusPV();

    virtual int allocate();
    virtual int tagDifferentiable();
    virtual int initialize();
    virtual int evaluateResidual();
    virtual int initializeAdjoint();
    virtual int evaluateAdjointResidual();

    virtual ScalarT& V()
    {
      return V_;
    }

    virtual const ScalarT& V() const
    {
      return V_;
    }

    virtual ScalarT& theta()
    {
      return y_[0];
    }

    virtual const ScalarT& theta() const
    {
      return y_[0];
    }

    virtual ScalarT& P()
    {
      return f_[0];
    }

    virtual const ScalarT& P() const
    {
      return f_[0];
    }

    virtual ScalarT& Q()
    {
      return Q_;
    }

    virtual const ScalarT& Q() const
    {
      return Q_;
    }

    virtual ScalarT& lambdaP()
    {
      assert(false);
      return thetaB_;
    }

    virtual const ScalarT& lambdaP() const
    {
      assert(false);
      return thetaB_;
    }

    virtual ScalarT& lambdaQ()
    {
      assert(false);
      return VB_;
    }

    virtual const ScalarT& lambdaQ() const
    {
      assert(false);
      return VB_;
    }

    virtual ScalarT& PB()
    {
      assert(false);
      return PB_;
    }

    virtual const ScalarT& PB() const
    {
      assert(false);
      return PB_;
    }

    virtual ScalarT& QB()
    {
      assert(false);
      return QB_;
    }

    virtual const ScalarT& QB() const
    {
      assert(false);
      return QB_;
    }

    virtual int BusType() const
    {
      return BaseBus<ScalarT, IdxT>::BusType::PV;
    }

  private:
    ScalarT V_;      ///< Bus voltage magnitude
    ScalarT theta0_; ///< Default initial value for phase
    ScalarT Q_;      ///< Reactive power that generator needs to provide

    ScalarT VB_;
    ScalarT thetaB_;
    ScalarT PB_;
    ScalarT QB_;
  };

} // namespace GridKit

#endif // _BUS_PV_HPP_
