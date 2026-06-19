
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
   * @brief Implementation of a fourth order generator model with
   * a simple governor.
   *
   */
  template <typename scalar_type, typename index_type>
  class Generator4Governor : public ModelEvaluatorImpl<scalar_type, index_type>
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

    Generator4Governor(BusT* bus, ScalarT P0, ScalarT Q0);
    virtual ~Generator4Governor();

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

  private:
    //
    // Private model methods
    //

    ScalarT Pg();
    ScalarT Qg();

    inline ScalarT frequencyPenalty(ScalarT omega);
    inline ScalarT frequencyPenaltyDer(ScalarT omega);

    inline ScalarT Lm(ScalarT Pm);
    inline ScalarT dLm(ScalarT Pm);

    inline ScalarT Ln(ScalarT Pn);
    inline ScalarT dLn(ScalarT Pn);

  public:
    //
    // Public inline accesor functions
    //

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

    const ScalarT& lambdaP() const
    {
      return bus_->lambdaP();
    }

    const ScalarT& lambdaQ() const
    {
      return bus_->lambdaQ();
    }

    ScalarT& PB()
    {
      return bus_->PB();
    }

    ScalarT& QB()
    {
      return bus_->QB();
    }

  private:
    //
    // Private inlined accessor methods
    //

    const ScalarT dotDelta() const
    {
      return yp_.getData()[static_cast<size_t>(offsetGen_ + 0)];
    }

    const ScalarT dotOmega() const
    {
      return yp_.getData()[static_cast<size_t>(offsetGen_ + 1)];
    }

    const ScalarT dotEdp() const
    {
      return yp_.getData()[static_cast<size_t>(offsetGen_ + 2)];
    }

    const ScalarT dotEqp() const
    {
      return yp_.getData()[static_cast<size_t>(offsetGen_ + 3)];
    }

    const ScalarT delta() const
    {
      return y_.getData()[static_cast<size_t>(offsetGen_ + 0)];
    }

    const ScalarT omega() const
    {
      return y_.getData()[static_cast<size_t>(offsetGen_ + 1)];
    }

    const ScalarT Edp() const
    {
      return y_.getData()[static_cast<size_t>(offsetGen_ + 2)];
    }

    const ScalarT Eqp() const
    {
      return y_.getData()[static_cast<size_t>(offsetGen_ + 3)];
    }

    const ScalarT Id() const
    {
      return y_.getData()[static_cast<size_t>(offsetGen_ + 4)];
    }

    const ScalarT Iq() const
    {
      return y_.getData()[static_cast<size_t>(offsetGen_ + 5)];
    }

    const ScalarT K() const
    {
      return param_.getData()[1];
    }

    const ScalarT T1() const
    {
      return T1_;
    }

    const ScalarT T2() const
    {
      return param_.getData()[0];
    }

    const ScalarT T3() const
    {
      return T3_;
    }

  private:
    // Generator parameters
    RealT H_;
    RealT D_;
    RealT Xq_;
    RealT Xd_;
    RealT Xqp_;
    RealT Xdp_;
    RealT Rs_;
    RealT Tq0p_;
    RealT Td0p_;
    RealT Ef0_;
    RealT Pm0_;
    RealT deltaPm_;
    RealT deltaPn_;
    RealT omega_s_;
    RealT omega_b_;
    RealT omega_up_;
    RealT omega_lo_;
    RealT c_;
    RealT beta_;

    // Governor parameters
    RealT T1_;
    RealT T2_;
    RealT T3_;
    RealT K_;

    // Index offsets
    const IdxT offsetGen_;
    const IdxT offsetGov_;

    // Initial power flow values
    ScalarT P0_;
    ScalarT Q0_;

    // Bus to which the generator is connected
    BusT* bus_;
  };

} // namespace GridKit
