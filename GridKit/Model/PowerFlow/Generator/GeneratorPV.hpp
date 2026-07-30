
#pragma once

#include <vector>

#include <GridKit/Model/PowerFlow/Generator/GeneratorBase.hpp>
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
   * @brief Implementation of a PV generator.
   *
   */
  template <typename scalar_type, typename index_type>
  class GeneratorPV : public GeneratorBase<scalar_type, index_type>
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

    GeneratorPV(BusT* bus, GenData& data);
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
    ScalarT P_;
    BusT*   bus_;
  };
} // namespace GridKit
