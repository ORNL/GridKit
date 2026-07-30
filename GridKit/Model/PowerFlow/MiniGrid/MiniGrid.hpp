#pragma once

#include <vector>

#include <GridKit/Model/PowerFlow/ModelEvaluatorImpl.hpp>

namespace GridKit
{
  /*!
   * @brief Implementation of a power grid.
   *
   */
  template <typename scalar_type, typename index_type>
  class MiniGrid : public ModelEvaluatorImpl<scalar_type, index_type>
  {
    using ModelEvaluatorImpl<scalar_type, index_type>::size_;
    using ModelEvaluatorImpl<scalar_type, index_type>::nnz_;
    using ModelEvaluatorImpl<scalar_type, index_type>::time_;
    using ModelEvaluatorImpl<scalar_type, index_type>::y_;
    using ModelEvaluatorImpl<scalar_type, index_type>::f_;

  public:
    using ScalarT = scalar_type;
    using IdxT    = index_type;
    using RealT   = typename ModelEvaluatorImpl<ScalarT, IdxT>::RealT;

    MiniGrid();
    virtual ~MiniGrid();

    int allocate();
    int initialize();

    int tagDifferentiable()
    {
      return -1;
    }

    int setAbsoluteTolerance(RealT)
    {
      return -1;
    }

    int evaluateResidual();
    int evaluateJacobian();

    int evaluateIntegrand()
    {
      return -1;
    }

    int initializeAdjoint()
    {
      return -1;
    }

    int evaluateAdjointResidual()
    {
      return -1;
    }

    // int evaluateAdjointJacobian() {return -1;}
    int evaluateAdjointIntegrand()
    {
      return -1;
    }

    void updateTime(RealT /* t */, RealT /* a */)
    {
    }

    // const accessors are public
    ScalarT const& th2() const
    {
      return y_.getData()[0];
    }

    ScalarT const& V2() const
    {
      return y_.getData()[1];
    }

    ScalarT const& th3() const
    {
      return y_.getData()[2];
    }

    ScalarT& th2()
    {
      return y_.getData()[0];
    }

    ScalarT& V2()
    {
      return y_.getData()[1];
    }

    ScalarT& th3()
    {
      return y_.getData()[2];
    }

  private:
    ScalarT Pl2_;
    ScalarT Ql2_;
    ScalarT Pg3_;
    ScalarT V1_;
    ScalarT th1_;
    ScalarT V3_;
    ScalarT B12_;
    ScalarT B13_;
    ScalarT B22_;
    ScalarT B23_;
  };
} // namespace GridKit
