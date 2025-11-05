
#pragma once

#include <GridKit/Model/PowerFlow/ModelEvaluatorImpl.hpp>

namespace GridKit
{
  /*!
   * @brief Base class for all power flow buses.
   *
   * Derived bus types:
   *   0 - swing bus (V and theta are constants)
   *   1 - PV bus    (P and V are constants)
   *   2 - PQ bus    (P and Q are constants)
   *
   * @todo Consider static instead of dynamic polymorphism for
   * bus types. Create Bus class that takes template parameter
   * BusType.
   */
  template <class ScalarT, typename IdxT>
  class BaseBus : public ModelEvaluatorImpl<ScalarT, IdxT>
  {
  protected:
    using ModelEvaluatorImpl<ScalarT, IdxT>::size_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::nnz_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::time_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::alpha_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::rel_tol_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::abs_tol_;
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

    enum BusType
    {
      PQ = 1,
      PV,
      Slack,
      Isolated
    };

    BaseBus(IdxT id)
      : busID_(id)
    {
    }

    virtual ~BaseBus()
    {
    }

    // Set defaults for ModelEvaluator methods
    virtual int allocate()
    {
      return 0;
    }

    virtual int initialize()
    {
      return 0;
    }

    virtual int tagDifferentiable()
    {
      return 0;
    }

    virtual int evaluateResidual()
    {
      return 0;
    }

    virtual int evaluateJacobian()
    {
      return 0;
    }

    virtual int evaluateIntegrand()
    {
      return 0;
    }

    virtual int initializeAdjoint()
    {
      return 0;
    }

    virtual int evaluateAdjointResidual()
    {
      return 0;
    }

    // virtual int evaluateAdjointJacobian() { return 0;}
    virtual int evaluateAdjointIntegrand()
    {
      return 0;
    }

    virtual void updateTime(real_type, real_type)
    {
    } // <- throw exception here

    // Pure virtual methods specific to Bus types
    virtual ScalarT&       V()           = 0;
    virtual const ScalarT& V() const     = 0;
    virtual ScalarT&       theta()       = 0;
    virtual const ScalarT& theta() const = 0;
    virtual ScalarT&       P()           = 0;
    virtual const ScalarT& P() const     = 0;
    virtual ScalarT&       Q()           = 0;
    virtual const ScalarT& Q() const     = 0;

    virtual ScalarT&       lambdaP()       = 0;
    virtual const ScalarT& lambdaP() const = 0;
    virtual ScalarT&       lambdaQ()       = 0;
    virtual const ScalarT& lambdaQ() const = 0;
    virtual ScalarT&       PB()            = 0;
    virtual const ScalarT& PB() const      = 0;
    virtual ScalarT&       QB()            = 0;
    virtual const ScalarT& QB() const      = 0;

    virtual int BusType() const = 0;

    virtual const IdxT BusID() const
    {
      return busID_;
    }

  protected:
    const IdxT busID_;
  }; // class BaseBus

} // namespace GridKit
