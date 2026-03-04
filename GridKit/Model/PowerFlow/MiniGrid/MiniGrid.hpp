#pragma once

#include <vector>

#include <GridKit/Model/PowerFlow/ModelEvaluatorImpl.hpp>

namespace GridKit
{
  /*!
   * @brief Implementation of a power grid.
   *
   */
  template <class ScalarT, typename IdxT>
  class MiniGrid : public ModelEvaluatorImpl<ScalarT, IdxT>
  {
    using ModelEvaluatorImpl<ScalarT, IdxT>::size_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::nnz_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::time_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::y_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::f_;

    using RealT = typename ModelEvaluatorImpl<ScalarT, IdxT>::RealT;

  public:
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
      return y_[0];
    }

    ScalarT const& V2() const
    {
      return y_[1];
    }

    ScalarT const& th3() const
    {
      return y_[2];
    }

    ScalarT& th2()
    {
      return y_[0];
    }

    ScalarT& V2()
    {
      return y_[1];
    }

    ScalarT& th3()
    {
      return y_[2];
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
