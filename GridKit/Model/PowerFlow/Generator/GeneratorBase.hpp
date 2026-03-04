
#pragma once

#include <vector>

#include <GridKit/Model/PowerFlow/ModelEvaluatorImpl.hpp>

namespace GridKit
{
  template <class ScalarT, typename IdxT>
  class BaseBus;
}

namespace GridKit
{
  /**
   * @brief Generator base class template
   *
   * @tparam ScalarT - Scalar type
   * @tparam IdxT    - Matrix and vector index type
   */
  template <class ScalarT, typename IdxT>
  class GeneratorBase : public ModelEvaluatorImpl<ScalarT, IdxT>
  {
  protected:
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

    using bus_type = BaseBus<ScalarT, IdxT>;
    using RealT    = typename ModelEvaluatorImpl<ScalarT, IdxT>::RealT;

  public:
    GeneratorBase()
    {
    }

    virtual ~GeneratorBase()
    {
    }

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

    // virtual int evaluateAdjointJacobian();
    virtual int evaluateAdjointIntegrand()
    {
      return 0;
    }

    void updateTime(RealT, RealT)
    {
    }

    virtual ScalarT&       P()       = 0;
    virtual const ScalarT& P() const = 0;
    virtual ScalarT&       Q()       = 0;
    virtual const ScalarT& Q() const = 0;
  };
} // namespace GridKit
