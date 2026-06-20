
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
    template <typename scalar_type, typename index_type>
    class BusInfinite : public BusBase<scalar_type, index_type>
    {
      using BusBase<scalar_type, index_type>::bus_id_;
      using BusBase<scalar_type, index_type>::size_;
      using BusBase<scalar_type, index_type>::y_;
      using BusBase<scalar_type, index_type>::yp_;
      using BusBase<scalar_type, index_type>::f_;
      using BusBase<scalar_type, index_type>::J_;
      using BusBase<scalar_type, index_type>::variable_indices_;
      using BusBase<scalar_type, index_type>::residual_indices_;
      using BusBase<scalar_type, index_type>::monitor_;

    public:
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename BusBase<ScalarT, IdxT>::RealT;
      using MonitorT   = typename BusBase<ScalarT, IdxT>::MonitorT;
      using ModelDataT = BusData<RealT, IdxT>;
      using BusTypeT   = typename BusData<RealT, IdxT>::BusType;

      BusInfinite();
      BusInfinite(ScalarT Vr, ScalarT Vi);
      BusInfinite(const ModelDataT& data);
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
