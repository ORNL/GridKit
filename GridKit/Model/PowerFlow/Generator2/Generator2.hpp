
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
    using ModelEvaluatorImpl<ScalarT, IdxT>::f_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::g_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::yB_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::ypB_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::fB_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::gB_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::param_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::param_up_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::param_lo_;

    typedef typename ModelEvaluatorImpl<ScalarT, IdxT>::real_type real_type;
    typedef BaseBus<ScalarT, IdxT>                                bus_type;

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

    void updateTime(real_type t, real_type a)
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
    real_type H_;
    real_type D_;
    real_type Pm_;
    real_type Xdp_;
    real_type Eqp_;
    real_type omega_s_;
    real_type omega_b_;
    real_type omega_up_;
    real_type omega_lo_;
    real_type theta_s_;
    real_type c_;
    real_type beta_;

    bus_type* bus_;
  };

} // namespace GridKit
