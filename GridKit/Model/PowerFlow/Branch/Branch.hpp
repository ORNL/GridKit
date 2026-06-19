
#pragma once

#include <GridKit/Model/PowerFlow/ModelEvaluatorImpl.hpp>

// Forward declarations.
namespace GridKit
{
  template <typename scalar_type, typename index_type>
  class BaseBus;

  namespace PowerFlowData
  {
    template <typename scalar_type, typename index_type>
    struct BranchData;
  }
} // namespace GridKit

namespace GridKit
{
  /*!
   * @brief Implementation of a pi-model branch between two buses.
   *
   */
  template <typename scalar_type, typename index_type>
  class Branch : public ModelEvaluatorImpl<scalar_type, index_type>
  {
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
    using ScalarT    = scalar_type;
    using IdxT       = index_type;
    using RealT      = typename ModelEvaluatorImpl<ScalarT, IdxT>::RealT;
    using BusT       = BaseBus<ScalarT, IdxT>;
    using BranchData = GridKit::PowerFlowData::BranchData<RealT, IdxT>;

    Branch(BusT* bus1, BusT* bus2);
    Branch(RealT R, RealT X, RealT G, RealT B, BusT* bus1, BusT* bus2);
    Branch(BusT* bus1, BusT* bus2, BranchData& data);
    virtual ~Branch();

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

    void updateTime(RealT /* t */, RealT /* a */)
    {
    }

  public:
    void setR(RealT R)
    {
      R_ = R;
    }

    void setX(RealT X)
    {
      // std::cout << "Setting X ...\n";
      X_ = X;
    }

    void setG(RealT G)
    {
      G_ = G;
    }

    void setB(RealT B)
    {
      B_ = B;
    }

  private:
    ScalarT& V1()
    {
      return bus1_->V();
    }

    ScalarT& theta1()
    {
      return bus1_->theta();
    }

    ScalarT& P1()
    {
      return bus1_->P();
    }

    ScalarT& Q1()
    {
      return bus1_->Q();
    }

    ScalarT& V2()
    {
      return bus2_->V();
    }

    ScalarT& theta2()
    {
      return bus2_->theta();
    }

    ScalarT& P2()
    {
      return bus2_->P();
    }

    ScalarT& Q2()
    {
      return bus2_->Q();
    }

  private:
    RealT      R_;
    RealT      X_;
    RealT      G_;
    RealT      B_;
    const IdxT fbusID_;
    const IdxT tbusID_;
    BusT*      bus1_;
    BusT*      bus2_;
  };
} // namespace GridKit
