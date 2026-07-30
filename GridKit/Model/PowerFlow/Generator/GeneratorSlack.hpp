
#pragma once

#include <vector>

#include <GridKit/Model/PowerFlow/Generator/GeneratorBase.hpp>
#include <GridKit/Model/PowerFlow/PowerFlowData.hpp>

namespace GridKit
{
  template <typename scalar_type, typename index_type>
  class BaseBus;
}

namespace GridKit
{
  /*!
   * @brief Implementation of a power grid.
   *
   */
  template <typename scalar_type, typename index_type>
  class GeneratorSlack : public GeneratorBase<scalar_type, index_type>
  {
    using GeneratorBase<scalar_type, index_type>::size_;
    using GeneratorBase<scalar_type, index_type>::nnz_;
    using GeneratorBase<scalar_type, index_type>::time_;
    using GeneratorBase<scalar_type, index_type>::alpha_;
    using GeneratorBase<scalar_type, index_type>::y_;
    using GeneratorBase<scalar_type, index_type>::yp_;
    using GeneratorBase<scalar_type, index_type>::tag_;
    using GeneratorBase<scalar_type, index_type>::f_;
    using GeneratorBase<scalar_type, index_type>::g_;
    using GeneratorBase<scalar_type, index_type>::yB_;
    using GeneratorBase<scalar_type, index_type>::ypB_;
    using GeneratorBase<scalar_type, index_type>::fB_;
    using GeneratorBase<scalar_type, index_type>::gB_;
    using GeneratorBase<scalar_type, index_type>::param_;

  public:
    using ScalarT = scalar_type;
    using IdxT    = index_type;
    using RealT   = typename ModelEvaluatorImpl<ScalarT, IdxT>::RealT;
    using BusT    = BaseBus<ScalarT, IdxT>;
    using GenData = GridKit::PowerFlowData::GenData<RealT, IdxT>;

    GeneratorSlack(BusT* bus, GenData& data);
    virtual ~GeneratorSlack();

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

    virtual ScalarT& P()
    {
      return bus_->P();
    }

    virtual const ScalarT& P() const
    {
      return bus_->P();
    }

    virtual ScalarT& Q()
    {
      return bus_->Q();
    }

    virtual const ScalarT& Q() const
    {
      return bus_->Q();
    }

  private:
    BusT* bus_;
  };
} // namespace GridKit
