
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
     * @brief Implementation of an "infinite" bus.
     *
     *
     *
     */
    template <class ScalarT, typename IdxT>
    class BusInfinite : public BusElectric<ScalarT, IdxT>
    {
      using BusElectric<ScalarT, IdxT>::size_;
      using BusElectric<ScalarT, IdxT>::y_;
      using BusElectric<ScalarT, IdxT>::yp_;
      using BusElectric<ScalarT, IdxT>::yB_;
      using BusElectric<ScalarT, IdxT>::ypB_;
      using BusElectric<ScalarT, IdxT>::f_;
      using BusElectric<ScalarT, IdxT>::fB_;
      using BusElectric<ScalarT, IdxT>::tag_;

    public:
      using real_type = typename BusElectric<ScalarT, IdxT>::real_type;
      using DataT     = BusData<real_type, IdxT>;

      BusInfinite();
      BusInfinite(ScalarT Vr, ScalarT Vi);
      BusInfinite(const DataT& data);
      virtual ~BusInfinite();

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
        return BusData<real_type, IdxT>::BusType::SLACK;
      }

    private:
      ScalarT Vr_{0.0};
      ScalarT Vi_{0.0};
      ScalarT Ir_{0.0};
      ScalarT Ii_{0.0};

      ScalarT VrB_{0.0};
      ScalarT ViB_{0.0};
      ScalarT IrB_{0.0};
      ScalarT IiB_{0.0};
    };

  } // namespace PhasorDynamics
} // namespace GridKit
