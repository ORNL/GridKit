
#pragma once

#include <vector>

#include <GridKit/Model/PowerFlow/Generator/GeneratorBase.hpp>
#include <GridKit/Model/PowerFlow/PowerFlowData.hpp>

namespace GridKit
{
  template <class ScalarT, typename IdxT>
  class BaseBus;
}

namespace GridKit
{
  /*!
   * @brief Implementation of a power grid.
   *
   */
  template <class ScalarT, typename IdxT>
  class GeneratorSlack : public GeneratorBase<ScalarT, IdxT>
  {
    using GeneratorBase<ScalarT, IdxT>::size_;
    using GeneratorBase<ScalarT, IdxT>::nnz_;
    using GeneratorBase<ScalarT, IdxT>::time_;
    using GeneratorBase<ScalarT, IdxT>::alpha_;
    using GeneratorBase<ScalarT, IdxT>::y_;
    using GeneratorBase<ScalarT, IdxT>::yp_;
    using GeneratorBase<ScalarT, IdxT>::tag_;
    using GeneratorBase<ScalarT, IdxT>::f_;
    using GeneratorBase<ScalarT, IdxT>::g_;
    using GeneratorBase<ScalarT, IdxT>::yB_;
    using GeneratorBase<ScalarT, IdxT>::ypB_;
    using GeneratorBase<ScalarT, IdxT>::fB_;
    using GeneratorBase<ScalarT, IdxT>::gB_;
    using GeneratorBase<ScalarT, IdxT>::param_;

    using bus_type = BaseBus<ScalarT, IdxT>;
    using RealT    = typename ModelEvaluatorImpl<ScalarT, IdxT>::RealT;
    using GenData  = GridKit::PowerFlowData::GenData<RealT, IdxT>;

  public:
    GeneratorSlack(bus_type* bus, GenData& data);
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
    bus_type* bus_;
  };
} // namespace GridKit
