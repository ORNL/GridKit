
#pragma once

#include <vector>

#include <GridKit/Model/PowerFlow/Generator/GeneratorBase.hpp>
#include <GridKit/Model/PowerFlow/ModelEvaluatorImpl.hpp>
#include <GridKit/Model/PowerFlow/PowerFlowData.hpp>

namespace GridKit
{
  template <class ScalarT, typename IdxT>
  class BaseBus;
}

namespace GridKit
{
  /*!
   * @brief Implementation of a PV generator.
   *
   */
  template <class ScalarT, typename IdxT>
  class GeneratorPV : public GeneratorBase<ScalarT, IdxT>
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
    GeneratorPV(bus_type* bus, GenData& data);
    virtual ~GeneratorPV();

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
      return P_;
    }

    virtual const ScalarT& P() const
    {
      return P_;
    }

    /// @brief  Reactive power excess on PV bus
    /// @return reference to negative PV generator reactive power
    virtual ScalarT& Q()
    {
      return (bus_->Q());
    }

    /// @brief  Reactive power excess on PV bus
    /// @return const reference to negative PV generator reactive power
    virtual const ScalarT& Q() const
    {
      return (bus_->Q());
    }

  private:
    ScalarT   P_;
    bus_type* bus_;
  };
} // namespace GridKit
