/**
 * @file BusSignalVoltageOut.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Declaration of a bus with signal ports.
 */

#pragma once

#include <GridKit/Constants.hpp>
#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/SignalPorts.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /*!
     * @brief Bus with signal ports.
     *
     * Like @ref Bus, this model owns the bus voltage components _Vr_ and
     * _Vi_ as algebraic variables and uses current balance in Cartesian
     * coordinates as residuals. In addition, it
     * - publishes _Vr_ and _Vi_ on signal outlets `vr` and `vi`, and
     * - sets its residuals f[0] and f[1] to the current injections read from
     *   signal inlets `ir` and `ii`, respectively. Both inlets are
     *   mandatory: verify() throws if either is not connected to a linked
     *   signal, and no default value is ever used.
     *
     * Components attached to the bus directly (without signals) keep adding
     * their currents to the residuals after the bus residual is evaluated,
     * exactly as they do for @ref Bus.
     *
     * @note Signal ports have to be connected before allocate() is called, since
     *       the output signals are linked to the bus variables there.
     */
    template <typename scalar_type, typename index_type>
    class BusSignalVoltageOut : public BusBase<scalar_type, index_type>
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
      using ScalarT      = scalar_type;
      using IdxT         = index_type;
      using RealT        = typename BusBase<ScalarT, IdxT>::RealT;
      using CooMatrixT   = typename BusBase<ScalarT, IdxT>::CooMatrixT;
      using MonitorT     = typename BusBase<ScalarT, IdxT>::MonitorT;
      using ModelDataT   = BusData<RealT, IdxT>;
      using BusTypeT     = typename BusData<RealT, IdxT>::BusType;
      using SignalPortsT = SignalPorts<ScalarT, ModelDataT>;

      BusSignalVoltageOut();
      BusSignalVoltageOut(ScalarT Vr, ScalarT Vi);
      BusSignalVoltageOut(const ModelDataT& data);
      virtual ~BusSignalVoltageOut();

      virtual int setBusID(IdxT) override final;
      virtual int allocate() override final;
      virtual int verify() const override final;
      virtual int tagDifferentiable() override final;
      virtual int setAbsoluteTolerance(RealT rel_tol) override final;
      virtual int initialize() override final;
      virtual int evaluateResidual() override final;
      virtual int evaluateJacobian() override final;

      virtual BusTypeT BusType() const override final
      {
        return BusTypeT::SIGNAL_VOLTAGE_OUT;
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

      void setVr(RealT vr) override final
      {
        vr_init_ = static_cast<ScalarT>(vr);
      }

      void setVi(RealT vi) override final
      {
        vi_init_ = static_cast<ScalarT>(vi);
      }

      SignalPortsT& getPorts()
      {
        return ports_;
      }

      const SignalPortsT& getPorts() const
      {
        return ports_;
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
          coo_jac_->setDataPointers(J_rows_buffer_, J_cols_buffer_, J_vals_buffer_, memory::HOST);
        }

        return 0;
      }

      /**
       * @brief Initialize DependencyTracking variable numbers.
       *
       * @note Assigns even indices to y and odd indices to yp.
       *       Should be called in initialize(), after variables have been set.
       */
      int initializeDependencyTrackingVariableNumbers()
        requires std::is_same_v<ScalarT, DependencyTracking::Variable>
      {
        auto* y  = y_.getData();
        auto* yp = yp_.getData();

        for (IdxT j = 0; j < size_; ++j)
        {
          const IdxT var_idx = this->getVariableIndex(j);
          if (var_idx != INVALID_INDEX<IdxT>)
          {
            // Even indices for y and odd indices for yp
            y[j].setVariableNumber(static_cast<size_t>(2 * var_idx));
            yp[j].setVariableNumber(static_cast<size_t>(2 * var_idx + 1));
          }
        }

        y_.setDataUpdated();
        yp_.setDataUpdated();

        return 0;
      }

      IdxT*  J_rows_buffer_{nullptr};
      IdxT*  J_cols_buffer_{nullptr};
      RealT* J_vals_buffer_{nullptr};

    private:
      ScalarT vr_init_{0.0};
      ScalarT vi_init_{0.0};

      /// Signal ports
      SignalPortsT ports_;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
