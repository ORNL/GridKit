
#pragma once

#include <Model/PhasorDynamics/Bus/BusElectric/BusElectric.hpp>

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
    class BusNetwork : public BusNetwork<ScalarT, IdxT>
    {
      using BusBase<ScalarT, IdxT>::size_;
      using BusBase<ScalarT, IdxT>::y_;
      using BusBase<ScalarT, IdxT>::yp_;
      using BusBase<ScalarT, IdxT>::yB_;
      using BusBase<ScalarT, IdxT>::ypB_;
      using BusBase<ScalarT, IdxT>::f_;
      using BusBase<ScalarT, IdxT>::fB_;
      using BusBase<ScalarT, IdxT>::tag_;

    public:
      using real_type = typename BusNetwork<ScalarT, IdxT>::real_type;
      using DataT     = BusData<real_type, IdxT>;

      BusNetwork();
      BusNetwork(ScalarT Vr, ScalarT Vi);
      BusNetwork(const DataT& data);
      virtual ~BusNetwork();

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
        return BusData<real_type, IdxT>::BusType::DEFAULT;
      }

    private:
      ScalarT Vr0_{0.0};
      ScalarT Vi0_{0.0};
    };

  } // namespace PhasorDynamics
} // namespace GridKit
