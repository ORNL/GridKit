
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
   * @brief Implementation of a fourth order generator model with
   * a simple governor.
   *
   */
  template <class ScalarT, typename IdxT>
  class Generator4Governor : public ModelEvaluatorImpl<ScalarT, IdxT>
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

  public:
    typedef typename ModelEvaluatorImpl<ScalarT, IdxT>::real_type real_type;
    typedef BaseBus<ScalarT, IdxT>                                bus_type;

    Generator4Governor(bus_type* bus, ScalarT P0, ScalarT Q0);
    virtual ~Generator4Governor();

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
      return yp_[static_cast<size_t>(offsetGen_ + 0)];
    }

    const ScalarT dotOmega() const
    {
      return yp_[static_cast<size_t>(offsetGen_ + 1)];
    }

    const ScalarT dotEdp() const
    {
      return yp_[static_cast<size_t>(offsetGen_ + 2)];
    }

    const ScalarT dotEqp() const
    {
      return yp_[static_cast<size_t>(offsetGen_ + 3)];
    }

    const ScalarT delta() const
    {
      return y_[static_cast<size_t>(offsetGen_ + 0)];
    }

    const ScalarT omega() const
    {
      return y_[static_cast<size_t>(offsetGen_ + 1)];
    }

    const ScalarT Edp() const
    {
      return y_[static_cast<size_t>(offsetGen_ + 2)];
    }

    const ScalarT Eqp() const
    {
      return y_[static_cast<size_t>(offsetGen_ + 3)];
    }

    const ScalarT Id() const
    {
      return y_[static_cast<size_t>(offsetGen_ + 4)];
    }

    const ScalarT Iq() const
    {
      return y_[static_cast<size_t>(offsetGen_ + 5)];
    }

    const ScalarT K() const
    {
      return param_[1];
    }

    const ScalarT T1() const
    {
      return T1_;
    }

    const ScalarT T2() const
    {
      return param_[0];
    }

    const ScalarT T3() const
    {
      return T3_;
    }

  private:
    // Generator parameters
    real_type H_;
    real_type D_;
    real_type Xq_;
    real_type Xd_;
    real_type Xqp_;
    real_type Xdp_;
    real_type Rs_;
    real_type Tq0p_;
    real_type Td0p_;
    real_type Ef0_;
    real_type Pm0_;
    real_type deltaPm_;
    real_type deltaPn_;
    real_type omega_s_;
    real_type omega_b_;
    real_type omega_up_;
    real_type omega_lo_;
    real_type c_;
    real_type beta_;

    // Governor parameters
    real_type T1_;
    real_type T2_;
    real_type T3_;
    real_type K_;

    // Index offsets
    const IdxT offsetGen_;
    const IdxT offsetGov_;

    // Initial power flow values
    ScalarT P0_;
    ScalarT Q0_;

    // Bus to which the generator is connected
    bus_type* bus_;
  };

} // namespace GridKit
