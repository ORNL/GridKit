
#pragma once

#include <Model/PhasorDynamics/Bus/BusElectric.hpp>


namespace GridKit
{
  namespace PhasorDynamics
  {
    /*!
     * @brief Implementation of an "infinite" bus.
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
      using DataT     = BusElectricData<real_type, IdxT>;

      BusInfinite();
      BusInfinite(ScalarT Vr, ScalarT Vi);
      BusInfinite(const DataT& data);
      virtual ~BusInfinite();

      int allocate() override;
      int tagDifferentiable() override;
      int initialize() override;
      int evaluateResidual() override;
      int evaluateIntegrand() override;
      int evaluateJacobian() override;
      int initializeAdjoint() override;
      int evaluateAdjointIntegrand() override;
      int evaluateAdjointResidual() override;

      int BusType() const override
      {
        return BusElectricData<real_type, IdxT>::BusType::SLACK;
      }
      
      virtual ScalarT& Vr()
      {
        return y_[0];
      }

      virtual const ScalarT& Vr() const
      {
        return y_[0];
      }

      virtual ScalarT& Vi()
      {
        return y_[1];
      }

      virtual const ScalarT& Vi() const
      {
        return y_[1];
      }

      virtual ScalarT& Ir()
      {
        return f_[0];
      }

      virtual const ScalarT& Ir() const
      {
        return f_[0];
      }

      virtual ScalarT& Ii()
      {
        return f_[1];
      }

      virtual const ScalarT& Ii() const
      {
        return f_[1];
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
