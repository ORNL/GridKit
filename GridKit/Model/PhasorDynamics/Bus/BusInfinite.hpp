
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
      virtual int setAbsoluteTolerance(RealT rel_tol) override final;
      virtual int initialize() override final;
      virtual int evaluateResidual() override final;
      virtual int evaluateJacobian() override final;

      virtual BusTypeT BusType() const override final
      {
        return BusTypeT::SLACK;
      }

      virtual ScalarT& Vr() override final
      {
        return Vr_;
      }

      virtual const ScalarT& Vr() const override final
      {
        return Vr_;
      }

      virtual ScalarT& Vi() override final
      {
        return Vi_;
      }

      virtual const ScalarT& Vi() const override final
      {
        return Vi_;
      }

      virtual ScalarT& Ir() override final
      {
        return Ir_;
      }

      virtual const ScalarT& Ir() const override final
      {
        return Ir_;
      }

      virtual ScalarT& Ii() override final
      {
        return Ii_;
      }

      virtual const ScalarT& Ii() const override final
      {
        return Ii_;
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
