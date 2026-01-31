
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
   * @brief Implementation of a fourth order generator model.
   *
   */
  template <class ScalarT, typename IdxT>
  class Generator4Param : public ModelEvaluatorImpl<ScalarT, IdxT>
  {
    using ModelEvaluatorImpl<ScalarT, IdxT>::size_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::nnz_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::time_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::alpha_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::y_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::yp_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::tag_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::absTol_;
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
    Generator4Param(BaseBus<ScalarT, IdxT>* bus, ScalarT P0 = 1.0, ScalarT Q0 = 0.0);
    virtual ~Generator4Param();

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

    // Inline accesor functions
    ScalarT& V()
    {
      return bus_->V();
    }

    const ScalarT& V() const
    {
      return bus_->V();
    }

    ScalarT& theta()
    {
      return bus_->theta();
    }

    const ScalarT& theta() const
    {
      return bus_->theta();
    }

    ScalarT& P()
    {
      return bus_->P();
    }

    const ScalarT& P() const
    {
      return bus_->P();
    }

    ScalarT& Q()
    {
      return bus_->Q();
    }

    const ScalarT& Q() const
    {
      return bus_->Q();
    }

    ScalarT trajectoryPenalty(ScalarT t) const;
    ScalarT trajectoryPenaltyDerEqp(ScalarT t) const;
    ScalarT trajectoryPenaltyDerEdp(ScalarT t) const;

    std::vector<std::vector<ScalarT>>& getLookupTable()
    {
      return table_;
    }

    std::vector<std::vector<ScalarT>> const& getLookupTable() const
    {
      return table_;
    }

  private:
    const ScalarT& H() const
    {
      return param_[0];
    }

    const ScalarT& Pm() const
    {
      return Pm_;
      // return param_[0];
    }

    const ScalarT& Ef() const
    {
      return Ef_;
      // return param_[1];
    }

    ScalarT Pg();
    ScalarT Qg();

  private:
    //
    // Private inlined accessor methods
    //

    const ScalarT dotDelta() const
    {
      return yp_[0];
    }

    const ScalarT dotOmega() const
    {
      return yp_[1];
    }

    const ScalarT dotEdp() const
    {
      return yp_[2];
    }

    const ScalarT dotEqp() const
    {
      return yp_[3];
    }

    const ScalarT delta() const
    {
      return y_[0];
    }

    const ScalarT omega() const
    {
      return y_[1];
    }

    const ScalarT Edp() const
    {
      return y_[2];
    }

    const ScalarT Eqp() const
    {
      return y_[3];
    }

    const ScalarT Id() const
    {
      return y_[4];
    }

    const ScalarT Iq() const
    {
      return y_[5];
    }

  private:
    RealT H_;    ///< Inertia constant [s]
    RealT D_;    ///< Damping constant [pu]
    RealT Xq_;   ///< q-axis synchronous reactance [pu]
    RealT Xd_;   ///< d-axis synchronous reactance [pu]
    RealT Xqp_;  ///< q-axis transient reactance [pu]
    RealT Xdp_;  ///< d-axis transient reactance [pu]
    RealT Rs_;   ///< stator armature resistance [pu]
    RealT Tq0p_; ///< q-axis open circuit transient time constant [s]
    RealT Td0p_; ///< d-axis open circuit transient time constant [s]
    RealT Ef_;
    RealT Pm_;
    RealT omega_s_;
    RealT omega_b_;

    ScalarT P0_;
    ScalarT Q0_;

    bus_type* bus_;

    /// Look-up table data. @todo This should be part of a separate model.
    std::vector<std::vector<ScalarT>> table_;
  };

} // namespace GridKit
