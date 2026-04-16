
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
      using BusBase<ScalarT, IdxT>::J_rows_buffer_;
      using BusBase<ScalarT, IdxT>::J_cols_buffer_;
      using BusBase<ScalarT, IdxT>::J_vals_buffer_;
      using BusBase<ScalarT, IdxT>::tag_;
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

      virtual int setBusID(IdxT) override final;
      virtual int allocate() override final;
      virtual int tagDifferentiable() override final;
      virtual int initialize() override final;
      virtual int evaluateResidual() override final;
      virtual int evaluateJacobian() override final;

      virtual BusTypeT BusType() const override final
      {
        return BusTypeT::DEFAULT;
      }

      virtual IdxT count() const override final
      {
        return 2;
      }

      virtual ScalarT& v(IdxT k) override final
      {
        return y_[static_cast<size_t>(k)];
      }

      virtual const ScalarT& v(IdxT k) const override final
      {
        return y_[static_cast<size_t>(k)];
      }

      virtual ScalarT& i(IdxT k) override final
      {
        return f_[static_cast<size_t>(k)];
      }

      virtual const ScalarT& i(IdxT k) const override final
      {
        return f_[static_cast<size_t>(k)];
      }

    private:
      ScalarT Vr0_{0.0};
      ScalarT Vi0_{0.0};
    };

  } // namespace PhasorDynamics
} // namespace GridKit
