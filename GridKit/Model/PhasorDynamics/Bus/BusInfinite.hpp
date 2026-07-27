
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
      using BusBase<scalar_type, index_type>::nnz_;
      using BusBase<scalar_type, index_type>::y_;
      using BusBase<scalar_type, index_type>::yp_;
      using BusBase<scalar_type, index_type>::f_;
      using BusBase<scalar_type, index_type>::variable_indices_;
      using BusBase<scalar_type, index_type>::residual_indices_;
      using BusBase<scalar_type, index_type>::coo_jac_;
      using BusBase<scalar_type, index_type>::monitor_;
      using BusBase<scalar_type, index_type>::allocated_;

    public:
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename BusBase<ScalarT, IdxT>::RealT;
      using CooMatrixT = typename BusBase<ScalarT, IdxT>::CooMatrixT;
      using MonitorT   = typename BusBase<ScalarT, IdxT>::MonitorT;
      using ModelDataT = BusData<RealT, IdxT>;
      using BusTypeT   = typename BusData<RealT, IdxT>::BusType;

      using BusBase<scalar_type, index_type>::Vr;
      using BusBase<scalar_type, index_type>::Vi;
      using BusBase<scalar_type, index_type>::Ir;
      using BusBase<scalar_type, index_type>::Ii;

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

    protected:
      int refreshTerminals() override final;

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
