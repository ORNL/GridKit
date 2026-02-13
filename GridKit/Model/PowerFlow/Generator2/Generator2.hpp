
#pragma once

#include <GridKit/Model/PowerFlow/ModelEvaluatorImpl.hpp>

namespace GridKit
{
  template <class ScalarT, typename IdxT>
  class BaseBus;
}

namespace GridKit
{
  /*!
   * @brief Implementation of a second order generator model.
   *
   */
  template <class ScalarT, typename IdxT>
  class Generator2 : public ModelEvaluatorImpl<ScalarT, IdxT>
  {
    using ModelEvaluatorImpl<ScalarT, IdxT>::size_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::nnz_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::time_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::alpha_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::y_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::yp_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::tag_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::abs_tol_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::f_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::g_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::yB_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::ypB_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::fB_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::gB_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::param_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::param_up_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::param_lo_;

    using RealT    = typename ModelEvaluatorImpl<ScalarT, IdxT>::RealT;
    using bus_type = BaseBus<ScalarT, IdxT>;

  public:
    Generator2(bus_type* bus);
    virtual ~Generator2();

    int allocate();
    int initialize();
    int tagDifferentiable();
    int evaluateResidual();
    int evaluateJacobian();
    int evaluateIntegrand();

    int initializeAdjoint();
    int evaluateAdjointResidual();
    // int evaluateAdjointJacobian();
    int evaluateAdjointIntegrand();

    void updateTime(RealT t, RealT a)
    {
      time_  = t;
      alpha_ = a;
    }

    const ScalarT& V() const
    {
      return bus_->V();
    }

    ScalarT& V()
    {
      return bus_->V();
    }

    const ScalarT& theta() const
    {
      return bus_->theta();
    }

    ScalarT& theta()
    {
      return bus_->theta();
    }

  private:
    inline ScalarT frequencyPenalty(ScalarT omega);
    inline ScalarT frequencyPenaltyDer(ScalarT omega);

  private:
    RealT H_;
    RealT D_;
    RealT Pm_;
    RealT Xdp_;
    RealT Eqp_;
    RealT omega_s_;
    RealT omega_b_;
    RealT omega_up_;
    RealT omega_lo_;
    RealT theta_s_;
    RealT c_;
    RealT beta_;

    bus_type* bus_;
  };

} // namespace GridKit
