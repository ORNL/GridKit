
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
  template <typename scalar_type, typename index_type>
  class BaseBus : public ModelEvaluatorImpl<scalar_type, index_type>
  {
  protected:
    using ModelEvaluatorImpl<scalar_type, index_type>::size_;
    using ModelEvaluatorImpl<scalar_type, index_type>::nnz_;
    using ModelEvaluatorImpl<scalar_type, index_type>::time_;
    using ModelEvaluatorImpl<scalar_type, index_type>::alpha_;
    using ModelEvaluatorImpl<scalar_type, index_type>::y_;
    using ModelEvaluatorImpl<scalar_type, index_type>::yp_;
    using ModelEvaluatorImpl<scalar_type, index_type>::tag_;
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

    virtual int setAbsoluteTolerance(RealT)
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

    virtual void updateTime(RealT, RealT)
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
