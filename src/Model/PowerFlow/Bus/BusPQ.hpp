
#ifndef _BUS_PQ_HPP_
#define _BUS_PQ_HPP_

#include "BaseBus.hpp"
#include <Model/PowerFlow/PowerSystemData.hpp>

namespace GridKit
{
  /*!
   * @brief Implementation of a PQ bus.
   *
   * Voltage _V_ and phase _theta_ are variables in PQ bus model.
   * Active and reactive power, _P_ and _Q_, are residual components.
   *
   *
   */
  template <class ScalarT, typename IdxT>
  class BusPQ : public BaseBus<ScalarT, IdxT>
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

    BusPQ();
    BusPQ(ScalarT V, ScalarT theta);
    BusPQ(BusData& data);
    virtual ~BusPQ();

    virtual int allocate();
    virtual int tagDifferentiable();
    virtual int initialize();
    virtual int evaluateResidual();
    virtual int initializeAdjoint();
    virtual int evaluateAdjointResidual();

    virtual ScalarT& V()
    {
      return y_[0];
    }

    virtual const ScalarT& V() const
    {
      return y_[0];
    }

    virtual ScalarT& theta()
    {
      return y_[1];
    }

    virtual const ScalarT& theta() const
    {
      return y_[1];
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
      return f_[1];
    }

    virtual const ScalarT& Q() const
    {
      return f_[1];
    }

    virtual ScalarT& lambdaP()
    {
      return yB_[0];
    }

    virtual const ScalarT& lambdaP() const
    {
      return yB_[0];
    }

    virtual ScalarT& lambdaQ()
    {
      return yB_[1];
    }

    virtual const ScalarT& lambdaQ() const
    {
      return yB_[1];
    }

    virtual ScalarT& PB()
    {
      return fB_[0];
    }

    virtual const ScalarT& PB() const
    {
      return fB_[0];
    }

    virtual ScalarT& QB()
    {
      return fB_[1];
    }

    virtual const ScalarT& QB() const
    {
      return fB_[1];
    }

    virtual int BusType() const
    {
      return BaseBus<ScalarT, IdxT>::BusType::PQ;
    }

  private:
    // Default initial values for voltage and phase on PQ bus
    ScalarT V0_;
    ScalarT theta0_;
  };

} // namespace GridKit

#endif // _BUS_PQ_HPP_
