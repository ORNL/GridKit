
#pragma once

#include <GridKit/Model/PowerFlow/ModelEvaluatorImpl.hpp>

namespace GridKit
{
  template <typename scalar_type, typename index_type>
  class BaseBus;
}

namespace GridKit
{
  /*!
   * @brief Implementation of a second order generator model.
   *
   */
  template <typename scalar_type, typename index_type>
  class Generator2 : public ModelEvaluatorImpl<scalar_type, index_type>
  {
    using ModelEvaluatorImpl<scalar_type, index_type>::size_;
    using ModelEvaluatorImpl<scalar_type, index_type>::nnz_;
    using ModelEvaluatorImpl<scalar_type, index_type>::time_;
    using ModelEvaluatorImpl<scalar_type, index_type>::alpha_;
    using ModelEvaluatorImpl<scalar_type, index_type>::y_;
    using ModelEvaluatorImpl<scalar_type, index_type>::yp_;
    using ModelEvaluatorImpl<scalar_type, index_type>::tag_;
    using ModelEvaluatorImpl<scalar_type, index_type>::abs_tol_;
    using ModelEvaluatorImpl<scalar_type, index_type>::f_;
    using ModelEvaluatorImpl<scalar_type, index_type>::g_;
    using ModelEvaluatorImpl<scalar_type, index_type>::yB_;
    using ModelEvaluatorImpl<scalar_type, index_type>::ypB_;
    using ModelEvaluatorImpl<scalar_type, index_type>::fB_;
    using ModelEvaluatorImpl<scalar_type, index_type>::gB_;
    using ModelEvaluatorImpl<scalar_type, index_type>::param_;
    using ModelEvaluatorImpl<scalar_type, index_type>::param_up_;
    using ModelEvaluatorImpl<scalar_type, index_type>::param_lo_;

  public:
    using ScalarT  = scalar_type;
    using IdxT     = index_type;
    using RealT    = typename ModelEvaluatorImpl<ScalarT, IdxT>::RealT;
    using BusT = BaseBus<ScalarT, IdxT>;

    Generator2(BusT* bus);
    virtual ~Generator2();

    int allocate();
    int initialize();
    int tagDifferentiable();
    int setAbsoluteTolerance(RealT);
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

    BusT* bus_;
  };

} // namespace GridKit
