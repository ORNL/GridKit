
#pragma once

#include <GridKit/Model/PhasorDynamics/BusBase.hpp>

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
    class Bus : public BusBase<ScalarT, IdxT>
    {
      using BusBase<ScalarT, IdxT>::bus_id_;
      using BusBase<ScalarT, IdxT>::size_;
      using BusBase<ScalarT, IdxT>::y_;
      using BusBase<ScalarT, IdxT>::yp_;
      using BusBase<ScalarT, IdxT>::f_;
      using BusBase<ScalarT, IdxT>::J_;
      using BusBase<ScalarT, IdxT>::tag_;
      using BusBase<ScalarT, IdxT>::absTol_;
      using BusBase<ScalarT, IdxT>::variable_indices_;
      using BusBase<ScalarT, IdxT>::residual_indices_;

    public:
      using RealT    = typename BusBase<ScalarT, IdxT>::RealT;
      using DataT    = BusData<RealT, IdxT>;
      using BusTypeT = typename BusData<RealT, IdxT>::BusType;

      Bus();
      Bus(ScalarT Vr, ScalarT Vi);
      Bus(const DataT& data);
      virtual ~Bus();

      virtual int setBusID(IdxT) override;
      virtual int allocate() override;
      virtual int tagDifferentiable() override;
      virtual int setAbsoluteTolerance() override;
      virtual int initialize() override;
      virtual int evaluateResidual() override;
      virtual int evaluateJacobian() override;

      virtual BusTypeT BusType() const override
      {
        return BusTypeT::DEFAULT;
      }

      virtual ScalarT& Vr() override
      {
        return y_[0];
      }

      virtual const ScalarT& Vr() const override
      {
        return y_[0];
      }

      virtual ScalarT& Vi() override
      {
        return y_[1];
      }

      virtual const ScalarT& Vi() const override
      {
        return y_[1];
      }

      virtual ScalarT& Ir() override
      {
        return f_[0];
      }

      virtual const ScalarT& Ir() const override
      {
        return f_[0];
      }

      virtual ScalarT& Ii() override
      {
        return f_[1];
      }

      virtual const ScalarT& Ii() const override
      {
        return f_[1];
      }

    private:
      ScalarT Vr0_{0.0};
      ScalarT Vi0_{0.0};
    };

  } // namespace PhasorDynamics
} // namespace GridKit
