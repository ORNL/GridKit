
#ifndef _BUS_SLACK_HPP_
#define _BUS_SLACK_HPP_

#include "BaseBus.hpp"
#include <PowerSystemData.hpp>

namespace GridKit
{
  /*!
   * @brief Implementation of a slack bus.
   *
   * Slack bus sets voltage _V_ and phase _theta_ as constants.
   * Active and reactive power, _P_ and _Q_, are component model outputs,
   * but are computed outside the BusSlack class.
   *
   *
   */
  template <class ScalarT, typename IdxT>
  class BusSlack : public BaseBus<ScalarT, IdxT>
  {
    using BaseBus<ScalarT, IdxT>::size_;
    using BaseBus<ScalarT, IdxT>::y_;
    using BaseBus<ScalarT, IdxT>::yp_;
    using BaseBus<ScalarT, IdxT>::f_;
    using BaseBus<ScalarT, IdxT>::g_;
    using BaseBus<ScalarT, IdxT>::abs_tol_;
    using BaseBus<ScalarT, IdxT>::rel_tol_;

  public:
    using real_type = typename ModelEvaluatorImpl<ScalarT, IdxT>::real_type;
    using BusData   = GridKit::PowerSystemData::BusData<real_type, IdxT>;

    BusSlack();
    BusSlack(ScalarT V, ScalarT theta);
    BusSlack(BusData& data);
    virtual ~BusSlack();
    virtual int evaluateResidual();
    virtual int evaluateAdjointResidual();

    /// @todo Should slack bus allow changing voltage?
    virtual ScalarT& V()
    {
      return V_;
    }

    virtual const ScalarT& V() const
    {
      return V_;
    }

    /// @todo Should slack bus allow changing phase?
    virtual ScalarT& theta()
    {
      return theta_;
    }

    virtual const ScalarT& theta() const
    {
      return theta_;
    }

    virtual ScalarT& P()
    {
      return P_;
    }

    virtual const ScalarT& P() const
    {
      return P_;
    }

    virtual ScalarT& Q()
    {
      return Q_;
    }

    virtual const ScalarT& Q() const
    {
      return Q_;
    }

    /// @todo Should slack bus allow changing voltage?
    virtual ScalarT& lambdaP()
    {
      return thetaB_;
    }

    virtual const ScalarT& lambdaP() const
    {
      return thetaB_;
    }

    /// @todo Should slack bus allow changing phase?
    virtual ScalarT& lambdaQ()
    {
      return VB_;
    }

    virtual const ScalarT& lambdaQ() const
    {
      return VB_;
    }

    virtual ScalarT& PB()
    {
      return PB_;
    }

    virtual const ScalarT& PB() const
    {
      return PB_;
    }

    virtual ScalarT& QB()
    {
      return QB_;
    }

    virtual const ScalarT& QB() const
    {
      return QB_;
    }

    virtual int BusType() const
    {
      return BaseBus<ScalarT, IdxT>::BusType::Slack;
    }

  private:
    ScalarT V_;
    ScalarT theta_;
    ScalarT P_;
    ScalarT Q_;

    ScalarT VB_;
    ScalarT thetaB_;
    ScalarT PB_;
    ScalarT QB_;

  }; // class BusSlack

} // namespace GridKit

#endif // _BUS_SLACK_HPP_
