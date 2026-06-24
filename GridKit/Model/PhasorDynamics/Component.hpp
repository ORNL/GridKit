#pragma once

#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Utilities/Errors.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    using Log = ::GridKit::Utilities::Logger;

    /**
     * @brief Component model implementation base class.
     */
    template <class scalar_type, typename index_type>
    class Component : public Model::Evaluator<scalar_type, index_type>
    {
    public:
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename Model::Evaluator<ScalarT, IdxT>::RealT;
      using CsrMatrixT = typename Model::Evaluator<ScalarT, IdxT>::CsrMatrixT;
      using CooMatrixT = typename Model::Evaluator<ScalarT, IdxT>::CooMatrixT;

      Component() = default;

      virtual ~Component()
      {
        if (J_rows_buffer_ != nullptr)
        {
          delete[] J_rows_buffer_;
          delete[] J_cols_buffer_;
          delete[] J_vals_buffer_;
          J_rows_buffer_ = nullptr;
          J_cols_buffer_ = nullptr;
          J_vals_buffer_ = nullptr;
        }

        if (coo_jac_ != nullptr)
        {
          delete coo_jac_;
          coo_jac_ = nullptr;
        }

        if (csr_jac_ != nullptr)
        {
          delete csr_jac_;
          csr_jac_ = nullptr;
        }

        if (map_to_csr_ != nullptr)
        {
          delete[] map_to_csr_;
          map_to_csr_ = nullptr;
        }
      }

      virtual int verify() const = 0;

      IdxT size() override final
      {
        return size_;
      }

      IdxT nnz() override final
      {
        return nnz_;
      }

      std::vector<ScalarT>& y() override
      {
        return y_;
      }

      const std::vector<ScalarT>& y() const override
      {
        return y_;
      }

      std::vector<ScalarT>& yp() override
      {
        return yp_;
      }

      const std::vector<ScalarT>& yp() const override
      {
        return yp_;
      }

      std::vector<bool>& tag() override
      {
        return tag_;
      }

      const std::vector<bool>& tag() const override
      {
        return tag_;
      }

      std::vector<ScalarT>& absoluteTolerance() override
      {
        return abs_tol_;
      }

      const std::vector<ScalarT>& absoluteTolerance() const override
      {
        return abs_tol_;
      }

      std::vector<ScalarT>& getResidual() override
      {
        return f_;
      }

      const std::vector<ScalarT>& getResidual() const override
      {
        return f_;
      }

      int setVariableIndex(IdxT local_index, IdxT global_index)
      {
        variable_indices_[static_cast<size_t>(local_index)] = global_index;
        return 0;
      }

      IdxT& getVariableIndex(IdxT local_index)
      {
        return variable_indices_[static_cast<size_t>(local_index)];
      }

      const std::vector<IdxT>& getVariableIndices() const
      {
        return variable_indices_;
      }

      int setResidualIndex(IdxT local_index, IdxT global_index)
      {
        residual_indices_[static_cast<size_t>(local_index)] = global_index;
        return 0;
      }

      IdxT& getResidualIndex(IdxT local_index)
      {
        return residual_indices_[static_cast<size_t>(local_index)];
      }

      const std::vector<IdxT>& getResidualIndices() const
      {
        return residual_indices_;
      }

      CsrMatrixT* getCsrJacobian() const override
      {
        return csr_jac_;
      }

      CooMatrixT* getCooJacobian() const
      {
        return coo_jac_;
      }

      /// @todo Remove this method. It should be part of DynamicSolver class.
      bool hasJacobian() override
      {
        return true;
      }

      void updateTime(RealT t, RealT a) override
      {
        time_  = t;
        alpha_ = a;
      }

      /**
       * @brief Set system frequency and power bases.
       *
       * @param[in] freq_system_base - System frequency base in Hz.
       * @param[in] va_system_base - System power base in VA.
       */
      void setSystemBase(RealT freq_system_base, RealT va_system_base)
      {
        freq_system_base_ = freq_system_base;
        va_system_base_   = va_system_base;
      }

      virtual int setGridKitComponentID(IdxT) = 0;

      IdxT getGridKitComponentID() const
      {
        return gridkit_component_id_;
      }

      int constructCsr()
      {
        if (coo_jac_ == nullptr)
        {
          constructCoo();
        }

        if (csr_jac_ == nullptr)
        {
          IdxT* row_ptrs = coo_jac_->getCsrRowData();

          nnz_ = coo_jac_->getNnz();

          IdxT*  cols = new IdxT[static_cast<size_t>(nnz_)];
          RealT* vals = new RealT[static_cast<size_t>(nnz_)];

          std::copy(coo_jac_->getColData(), coo_jac_->getColData() + nnz_, cols);
          std::copy(coo_jac_->getValues(), coo_jac_->getValues() + nnz_, vals);

          csr_jac_ = new CsrMatrixT(coo_jac_->getNumRows(), coo_jac_->getNumColumns(), nnz_, &row_ptrs, &cols, &vals);
        }

        return 0;
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

      IdxT              size_{0};
      IdxT              nnz_{0};
      /// Global (system-level) variable indices
      std::vector<IdxT> variable_indices_;
      /// Global (system-level) residual indices
      std::vector<IdxT> residual_indices_;

      std::vector<ScalarT> y_;
      std::vector<ScalarT> yp_;
      std::vector<bool>    tag_;
      std::vector<ScalarT> abs_tol_;
      std::vector<ScalarT> f_;
      std::vector<ScalarT> g_;

      IdxT*       J_rows_buffer_{nullptr};
      IdxT*       J_cols_buffer_{nullptr};
      RealT*      J_vals_buffer_{nullptr};
      IdxT*       map_to_csr_{nullptr};
      CsrMatrixT* csr_jac_{nullptr};
      CooMatrixT* coo_jac_{nullptr};

      //
      // Adjoint sensitivity members
      //

      std::vector<ScalarT> yB_{};
      std::vector<ScalarT> ypB_{};
      std::vector<ScalarT> fB_{};
      std::vector<ScalarT> gB_{};

      std::vector<ScalarT> param_{};
      std::vector<ScalarT> param_up_{};
      std::vector<ScalarT> param_lo_{};

      IdxT gridkit_component_id_{0};

      std::vector<ScalarT> wb_;
      std::vector<ScalarT> h_;

      RealT time_;
      RealT alpha_;

      RealT freq_system_base_{60.0};
      RealT va_system_base_{100.0e6};

      using NotImplementedError = GridKit::Utilities::NotImplementedError;

    public:
      // TODO: evaluate how this complies with xSDK guidelines

      [[noreturn]] IdxT sizeQuadrature() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] IdxT sizeParams() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] std::vector<ScalarT>& yB() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const std::vector<ScalarT>& yB() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] std::vector<ScalarT>& ypB() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const std::vector<ScalarT>& ypB() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] std::vector<ScalarT>& param() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const std::vector<ScalarT>& param() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] std::vector<ScalarT>& param_up() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const std::vector<ScalarT>& param_up() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] std::vector<ScalarT>& param_lo() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const std::vector<ScalarT>& param_lo() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] int evaluateIntegrand() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] int initializeAdjoint() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] int evaluateAdjointResidual() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] int evaluateAdjointIntegrand() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] std::vector<ScalarT>& getIntegrand() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const std::vector<ScalarT>& getIntegrand() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] std::vector<ScalarT>& getAdjointResidual() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const std::vector<ScalarT>& getAdjointResidual() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] std::vector<ScalarT>& getAdjointIntegrand() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const std::vector<ScalarT>& getAdjointIntegrand() const override
      {
        throw NotImplementedError(__func__);
      }
    };

  } // namespace PhasorDynamics
} // namespace GridKit
