
#pragma once

#include <GridKit/Model/PowerFlow/Bus/BaseBus.hpp>
#include <GridKit/Model/PowerFlow/PowerFlowData.hpp>

namespace GridKit
{
  /*!
   * @brief Implementation of a PQ bus.
   *
   * Voltage _V_ and phase _theta_ are variables in PQ bus model.
   * Active and reactive power, _P_ and _Q_, are residual components.
   *
   *
   */
  template <typename scalar_type, typename index_type>
  class BusPQ : public BaseBus<scalar_type, index_type>
  {
    using BaseBus<scalar_type, index_type>::size_;
    using BaseBus<scalar_type, index_type>::y_;
    using BaseBus<scalar_type, index_type>::yp_;
    using BaseBus<scalar_type, index_type>::yB_;
    using BaseBus<scalar_type, index_type>::ypB_;
    using BaseBus<scalar_type, index_type>::f_;
    using BaseBus<scalar_type, index_type>::fB_;
    using BaseBus<scalar_type, index_type>::tag_;
    using BaseBus<scalar_type, index_type>::abs_tol_;

  public:
    using ScalarT = scalar_type;
    using IdxT    = index_type;
    using RealT   = typename ModelEvaluatorImpl<ScalarT, IdxT>::RealT;
    using BusData = GridKit::PowerFlowData::BusData<RealT, IdxT>;

    BusPQ();
    BusPQ(ScalarT V, ScalarT theta);
    BusPQ(BusData& data);
    virtual ~BusPQ();

    virtual int allocate();
    virtual int tagDifferentiable();
    virtual int setAbsoluteTolerance(RealT);
    virtual int initialize();
    virtual int evaluateResidual();
    virtual int initializeAdjoint();
    virtual int evaluateAdjointResidual();

    virtual ScalarT& V()
    {
      return y_.getData()[0];
    }

    virtual const ScalarT& V() const
    {
      return y_.getData()[0];
    }

    virtual ScalarT& theta()
    {
      return y_.getData()[1];
    }

    virtual const ScalarT& theta() const
    {
      return y_.getData()[1];
    }

    virtual ScalarT& P()
    {
      return f_.getData()[0];
    }

    virtual const ScalarT& P() const
    {
      return f_.getData()[0];
    }

    virtual ScalarT& Q()
    {
      return f_.getData()[1];
    }

    virtual const ScalarT& Q() const
    {
      return f_.getData()[1];
    }

    virtual ScalarT& lambdaP()
    {
      return yB_.getData()[0];
    }

    virtual const ScalarT& lambdaP() const
    {
      return yB_.getData()[0];
    }

    virtual ScalarT& lambdaQ()
    {
      return yB_.getData()[1];
    }

    virtual const ScalarT& lambdaQ() const
    {
      return yB_.getData()[1];
    }

    virtual ScalarT& PB()
    {
      return fB_.getData()[0];
    }

    virtual const ScalarT& PB() const
    {
      return fB_.getData()[0];
    }

    virtual ScalarT& QB()
    {
      return fB_.getData()[1];
    }

    virtual const ScalarT& QB() const
    {
      return fB_.getData()[1];
    }

    virtual int BusType() const
    {
      return BaseBus<ScalarT, IdxT>::BusType::PQ;
    }

  private:
    // Default initial values for voltage and phase on PQ bus
    ScalarT V0_;
    ScalarT theta0_;
  };

} // namespace GridKit
