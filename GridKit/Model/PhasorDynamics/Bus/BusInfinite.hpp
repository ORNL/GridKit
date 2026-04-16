
#pragma once

#include <GridKit/Model/PhasorDynamics/BusBase.hpp>

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
    class BusInfinite : public BusBase<ScalarT, IdxT>
    {
      using BusBase<ScalarT, IdxT>::bus_id_;
      using BusBase<ScalarT, IdxT>::size_;
      using BusBase<ScalarT, IdxT>::y_;
      using BusBase<ScalarT, IdxT>::yp_;
      using BusBase<ScalarT, IdxT>::f_;
      using BusBase<ScalarT, IdxT>::J_;
      using BusBase<ScalarT, IdxT>::variable_indices_;
      using BusBase<ScalarT, IdxT>::residual_indices_;

    public:
      using RealT    = typename BusBase<ScalarT, IdxT>::RealT;
      using DataT    = BusData<RealT, IdxT>;
      using BusTypeT = typename BusData<RealT, IdxT>::BusType;

      BusInfinite();
      BusInfinite(ScalarT Vr, ScalarT Vi);
      BusInfinite(const DataT& data);
      virtual ~BusInfinite();

      virtual int setBusID(IdxT) override final;
      virtual int allocate() override final;
      virtual int tagDifferentiable() override final;
      virtual int initialize() override final;
      virtual int evaluateResidual() override final;
      virtual int evaluateJacobian() override final;

      virtual BusTypeT BusType() const override final
      {
        return BusTypeT::SLACK;
      }

      virtual IdxT count() const override final
      {
        return 2;
      }

      virtual ScalarT& v(IdxT k) override final
      {
        return v_[static_cast<size_t>(k)];
      }

      virtual const ScalarT& v(IdxT k) const override final
      {
        return v_[static_cast<size_t>(k)];
      }

      virtual ScalarT& i(IdxT k) override final
      {
        return i_[static_cast<size_t>(k)];
      }

      virtual const ScalarT& i(IdxT k) const override final
      {
        return i_[static_cast<size_t>(k)];
      }

    private:
      ScalarT v_[2]{};
      ScalarT i_[2]{};
    };

  } // namespace PhasorDynamics
} // namespace GridKit
