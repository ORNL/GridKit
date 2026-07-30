
#pragma once

#include <cassert>

#include <GridKit/Model/PowerFlow/Bus/BaseBus.hpp>
#include <GridKit/Model/PowerFlow/PowerFlowData.hpp>

namespace GridKit
{
  /*!
   * @brief Implementation of a PV bus.
   *
   * Voltage _V_ and phase _theta_ are variables in PV bus model.
   * Active and reactive power, _P_ and _Q_, are residual components.
   *
   *
   */
  template <typename scalar_type, typename index_type>
  class BusPV : public BaseBus<scalar_type, index_type>
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

    BusPV();
    BusPV(ScalarT V, ScalarT theta0);
    BusPV(BusData& data);
    virtual ~BusPV();

    virtual int allocate();
    virtual int tagDifferentiable();
    virtual int setAbsoluteTolerance(RealT);
    virtual int initialize();
    virtual int evaluateResidual();
    virtual int initializeAdjoint();
    virtual int evaluateAdjointResidual();

    virtual ScalarT& V()
    {
      return V_;
    }

    virtual const ScalarT& V() const
    {
      return V_;
    }

    virtual ScalarT& theta()
    {
      return y_.getData()[0];
    }

    virtual const ScalarT& theta() const
    {
      return y_.getData()[0];
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
      return Q_;
    }

    virtual const ScalarT& Q() const
    {
      return Q_;
    }

    virtual ScalarT& lambdaP()
    {
      assert(false);
      return thetaB_;
    }

    virtual const ScalarT& lambdaP() const
    {
      assert(false);
      return thetaB_;
    }

    virtual ScalarT& lambdaQ()
    {
      assert(false);
      return VB_;
    }

    virtual const ScalarT& lambdaQ() const
    {
      assert(false);
      return VB_;
    }

    virtual ScalarT& PB()
    {
      assert(false);
      return PB_;
    }

    virtual const ScalarT& PB() const
    {
      assert(false);
      return PB_;
    }

    virtual ScalarT& QB()
    {
      assert(false);
      return QB_;
    }

    virtual const ScalarT& QB() const
    {
      assert(false);
      return QB_;
    }

    virtual int BusType() const
    {
      return BaseBus<ScalarT, IdxT>::BusType::PV;
    }

  private:
    ScalarT V_;      ///< Bus voltage magnitude
    ScalarT theta0_; ///< Default initial value for phase
    ScalarT Q_;      ///< Reactive power that generator needs to provide

    ScalarT VB_;
    ScalarT thetaB_;
    ScalarT PB_;
    ScalarT QB_;
  };

} // namespace GridKit
