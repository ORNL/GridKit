
#pragma once

#include <GridKit/Model/PowerFlow/ModelEvaluatorImpl.hpp>
#include <GridKit/Model/PowerFlow/PowerFlowData.hpp>

namespace GridKit
{
  template <typename scalar_type, typename index_type>
  class BaseBus;
}

namespace GridKit
{
  /*!
   * @brief Declaration of a passive load class.
   *
   */
  template <typename scalar_type, typename index_type>
  class Load : public ModelEvaluatorImpl<scalar_type, index_type>
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
    using ScalarT  = scalar_type;
    using IdxT     = index_type;
    using RealT    = typename ModelEvaluatorImpl<ScalarT, IdxT>::RealT;
    using BusT     = BaseBus<ScalarT, IdxT>;
    using LoadData = GridKit::PowerFlowData::LoadData<RealT, IdxT>;

    Load(BusT* bus, ScalarT P, ScalarT Q);
    Load(BusT* bus, LoadData& data);
    virtual ~Load();

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
    ScalarT    P_;
    ScalarT    Q_;
    const IdxT busID_;
    BusT*      bus_;
  };
} // namespace GridKit
