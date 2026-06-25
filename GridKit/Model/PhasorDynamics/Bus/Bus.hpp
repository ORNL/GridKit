
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
    template <typename scalar_type, typename index_type>
    class Bus : public BusBase<scalar_type, index_type>
    {
      using BusBase<scalar_type, index_type>::bus_id_;
      using BusBase<scalar_type, index_type>::size_;
      using BusBase<scalar_type, index_type>::nnz_;
      using BusBase<scalar_type, index_type>::y_;
      using BusBase<scalar_type, index_type>::yp_;
      using BusBase<scalar_type, index_type>::f_;
      using BusBase<scalar_type, index_type>::tag_;
      using BusBase<scalar_type, index_type>::abs_tol_;
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

      Bus();
      Bus(ScalarT Vr, ScalarT Vi);
      Bus(const ModelDataT& data);
      virtual ~Bus();

      virtual int setBusID(IdxT) override final;
      virtual int allocate() override final;
      virtual int tagDifferentiable() override final;
      virtual int setAbsoluteTolerance(RealT rel_tol) override final;
      virtual int initialize() override final;
      virtual int evaluateResidual() override final;
      virtual int evaluateJacobian() override final;

      virtual BusTypeT BusType() const override final
      {
        return BusTypeT::DEFAULT;
      }

      virtual ScalarT& Vr() override final
      {
        return y_.getData()[0];
      }

      virtual const ScalarT& Vr() const override final
      {
        return y_.getData()[0];
      }

      virtual ScalarT& Vi() override final
      {
        return y_.getData()[1];
      }

      virtual const ScalarT& Vi() const override final
      {
        return y_.getData()[1];
      }

      virtual ScalarT& Ir() override final
      {
        return f_.getData()[0];
      }

      virtual const ScalarT& Ir() const override final
      {
        return f_.getData()[0];
      }

      virtual ScalarT& Ii() override final
      {
        return f_.getData()[1];
      }

      virtual const ScalarT& Ii() const override final
      {
        return f_.getData()[1];
      }

    protected:
      int constructCoo()
      {
        if (coo_jac_ == nullptr)
        {
          IdxT num_rows = 0;
          IdxT num_cols = 0;
          for (IdxT i = 0; i < nnz_; ++i)
          {
            if (J_rows_buffer_[i] + 1 > num_rows)
            {
              num_rows = J_rows_buffer_[i] + 1;
            }
            if (J_cols_buffer_[i] + 1 > num_cols)
            {
              num_cols = J_cols_buffer_[i] + 1;
            }
          }
          coo_jac_ = new CooMatrixT(num_rows, num_cols, nnz_);
          coo_jac_->setDataPointers(J_rows_buffer_, J_cols_buffer_, J_vals_buffer_, LinearAlgebra::memory::HOST);
        }

        return 0;
      }

      IdxT*  J_rows_buffer_{nullptr};
      IdxT*  J_cols_buffer_{nullptr};
      RealT* J_vals_buffer_{nullptr};

    private:
      ScalarT Vr0_{0.0};
      ScalarT Vi0_{0.0};
    };

  } // namespace PhasorDynamics
} // namespace GridKit
