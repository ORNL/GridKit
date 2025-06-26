/**
 * @file BusSignal.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of BusSignal class.
 *
 */

#pragma once

#include <Model/PhasorDynamics/Bus/BusControl/BusControl.hpp>

// Forward declaration of BusData structure
namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename RealT, typename IdxT>
    struct BusData;
  }
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    /*!
     * @brief Implementation of a PQ bus.
     *
     * Voltage _V_ and phase _theta_ are variables in PQ bus model.
     * Active and reactive power, _P_ and _Q_, are residual components.
     *
     *
     */
    template <class ScalarT, typename IdxT>
    class BusSignal : public BusControl<ScalarT, IdxT>
    {
      using BusControl<ScalarT, IdxT>::size_;
      using BusControl<ScalarT, IdxT>::y_;
      using BusControl<ScalarT, IdxT>::yp_;
      using BusControl<ScalarT, IdxT>::yB_;
      using BusControl<ScalarT, IdxT>::ypB_;
      using BusControl<ScalarT, IdxT>::f_;
      using BusControl<ScalarT, IdxT>::fB_;
      using BusControl<ScalarT, IdxT>::tag_;

    public:
      using real_type = typename BusControl<ScalarT, IdxT>::real_type;
      using DataT     = BusData<real_type, IdxT>;

      BusSignal();
      BusSignal(const DataT& data);
      virtual ~BusSignal();

      virtual int allocate() override;
      virtual int tagDifferentiable() override;
      virtual int initialize() override;
      virtual int evaluateResidual() override;
      virtual int evaluateIntegrand() override;
      virtual int evaluateJacobian() override;
      virtual int initializeAdjoint() override;
      virtual int evaluateAdjointIntegrand() override;
      virtual int evaluateAdjointResidual() override;

      virtual int BusType() const override
      {
        return BusData<real_type, IdxT>::BusType::SIGNAL;
      }



    private:
      ScalarT channel;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
