
#pragma once

#include <GridKit/Model/PowerFlow/ModelEvaluatorImpl.hpp>

// Forward declarations.
namespace GridKit
{
  template <class ScalarT, typename IdxT>
  class BaseBus;

  namespace PowerFlowData
  {
    template <class ScalarT, typename IdxT>
    struct BranchData;
  }
} // namespace GridKit

namespace GridKit
{
  /*!
   * @brief Implementation of a pi-model branch between two buses.
   *
   */
  template <class ScalarT, typename IdxT>
  class Branch : public ModelEvaluatorImpl<ScalarT, IdxT>
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

    using bus_type   = BaseBus<ScalarT, IdxT>;
    using RealT      = typename ModelEvaluatorImpl<ScalarT, IdxT>::RealT;
    using BranchData = GridKit::PowerFlowData::BranchData<RealT, IdxT>;

  public:
    Branch(bus_type* bus1, bus_type* bus2);
    Branch(RealT R, RealT X, RealT G, RealT B, bus_type* bus1, bus_type* bus2);
    Branch(bus_type* bus1, bus_type* bus2, BranchData& data);
    virtual ~Branch();

    int allocate();
    int initialize();
    int tagDifferentiable();
    int setAbsoluteTolerance();
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
    bus_type*  bus1_;
    bus_type*  bus2_;
  };
} // namespace GridKit
