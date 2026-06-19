
#pragma once

#include <vector>

#include <GridKit/Model/PowerFlow/ModelEvaluatorImpl.hpp>

namespace GridKit
{
  template <typename scalar_type, typename index_type>
  class BaseBus;
}

namespace GridKit
{
  /**
   * @brief Generator base class template
   *
   * @tparam scalar_type - Scalar type
   * @tparam index_type    - Matrix and vector index type
   */
  template <typename scalar_type, typename index_type>
  class GeneratorBase : public ModelEvaluatorImpl<scalar_type, index_type>
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

  public:
    using ScalarT = scalar_type;
    using IdxT    = index_type;
    using RealT   = typename ModelEvaluatorImpl<ScalarT, IdxT>::RealT;
    using BusT    = BaseBus<ScalarT, IdxT>;

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
