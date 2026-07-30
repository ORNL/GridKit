
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
   * @brief Implementation of a fourth order generator model.
   *
   */
  template <typename scalar_type, typename index_type>
  class Generator4Param : public ModelEvaluatorImpl<scalar_type, index_type>
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
    using ScalarT = scalar_type;
    using IdxT    = index_type;
    using RealT   = typename ModelEvaluatorImpl<ScalarT, IdxT>::RealT;
    using BusT    = BaseBus<ScalarT, IdxT>;

    Generator4Param(BusT* bus, ScalarT P0 = 1.0, ScalarT Q0 = 0.0);
    virtual ~Generator4Param();

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
      return param_.getData()[0];
    }

    const ScalarT& Pm() const
    {
      return Pm_;
    }

    const ScalarT& Ef() const
    {
      return Ef_;
    }

    ScalarT Pg();
    ScalarT Qg();

  private:
    //
    // Private inlined accessor methods
    //

    const ScalarT dotDelta() const
    {
      return yp_.getData()[0];
    }

    const ScalarT dotOmega() const
    {
      return yp_.getData()[1];
    }

    const ScalarT dotEdp() const
    {
      return yp_.getData()[2];
    }

    const ScalarT dotEqp() const
    {
      return yp_.getData()[3];
    }

    const ScalarT delta() const
    {
      return y_.getData()[0];
    }

    const ScalarT omega() const
    {
      return y_.getData()[1];
    }

    const ScalarT Edp() const
    {
      return y_.getData()[2];
    }

    const ScalarT Eqp() const
    {
      return y_.getData()[3];
    }

    const ScalarT Id() const
    {
      return y_.getData()[4];
    }

    const ScalarT Iq() const
    {
      return y_.getData()[5];
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

    BusT* bus_;

    /// Look-up table data. @todo This should be part of a separate model.
    std::vector<std::vector<ScalarT>> table_;
  };

} // namespace GridKit
