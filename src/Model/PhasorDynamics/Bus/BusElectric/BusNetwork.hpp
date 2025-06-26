
#pragma once

#include <Model/PhasorDynamics/Bus/BusElectric.hpp>

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
    class BusNetwork : public BusElectric<ScalarT, IdxT>
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
      
      // Constructors
      BusNetwork();
      BusNetwork(ScalarT Vr, ScalarT Vi);
      BusNetwork(const DataT& data);

      // Default Destructor
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

      int BusType() const override
      {
        return BusElectricData<real_type, IdxT>::BusType::DEFAULT;
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
      ScalarT Vr0_{0.0};
      ScalarT Vi0_{0.0};
    };

  } // namespace PhasorDynamics
} // namespace GridKit
