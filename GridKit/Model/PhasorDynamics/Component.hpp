#pragma once

#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/CommonMath.hpp>
#include <GridKit/Constants.hpp>
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
      using VectorT    = typename Model::Evaluator<ScalarT, IdxT>::VectorT;

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

      // @todo This is much clearer if we write it this way
      virtual int evaluateInternalResidual()
      {
        return this->evaluateResidual();
      }

      virtual int evaluateExternalResidual()
      {
        return 0;
      }

      IdxT size() override final
      {
        return size_;
      }

      IdxT nnz() override final
      {
        return nnz_;
      }

      VectorT& y() override
      {
        return y_;
      }

      const VectorT& y() const override
      {
        return y_;
      }

      VectorT& yp() override
      {
        return yp_;
      }

      const VectorT& yp() const override
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

      VectorT& absoluteTolerance() override
      {
        return abs_tol_;
      }

      const VectorT& absoluteTolerance() const override
      {
        return abs_tol_;
      }

      VectorT& getResidual() override
      {
        return f_;
      }

      const VectorT& getResidual() const override
      {
        return f_;
      }

      /**
       * @brief Bind this component's state and residual vectors to the slice
       * [offset, offset + size()) of the system vectors.
       *
       * After binding, the component reads and writes system storage directly
       * and allocate() will not allocate local vector data. Rebinding is
       * allowed and refreshes the aliases, e.g. after system storage is
       * reallocated when the topology changes.
       *
       * Only HOST data is bound because PhasorDynamics currently evaluates
       * models on the CPU. Supporting DEVICE data would also require sharing
       * the matching DEVICE pointer and keeping the HOST and DEVICE copies in
       * sync. This bind operation does neither, so DEVICE access is unsupported.
       *
       * @param[in] y       - System state vector
       * @param[in] yp      - System state derivative vector
       * @param[in] f       - System residual vector
       * @param[in] abs_tol - System absolute tolerance vector
       * @param[in] offset  - Position of this component's slice in the system vectors
       *
       * @pre System vectors hold current HOST data of at least offset + size()
       * elements. This component's vectors are unallocated or already bound.
       * @post allocated_ is true and y_, yp_, f_, abs_tol_ alias system storage.
       *
       * @return 0 if successful, non-zero otherwise.
       */
      int bind(VectorT& y, VectorT& yp, VectorT& f, VectorT& abs_tol, IdxT offset)
      {
        if (y.getSize() < offset + size_
            || yp.getSize() < offset + size_
            || f.getSize() < offset + size_
            || abs_tol.getSize() < offset + size_)
        {
          Log::error() << "Component::bind - system vectors are smaller than "
                       << "offset + size = " << offset + size_ << "\n";
          return 1;
        }

        auto* y_data       = y.getData(memory::HOST);
        auto* yp_data      = yp.getData(memory::HOST);
        auto* f_data       = f.getData(memory::HOST);
        auto* abs_tol_data = abs_tol.getData(memory::HOST);

        if (y_data == nullptr || yp_data == nullptr
            || f_data == nullptr || abs_tol_data == nullptr)
        {
          Log::error() << "Component::bind - system vector data is null or stale\n";
          return 1;
        }

        const int y_status       = y_.setData(y_data + offset, size_, memory::HOST);
        const int yp_status      = yp_.setData(yp_data + offset, size_, memory::HOST);
        const int f_status       = f_.setData(f_data + offset, size_, memory::HOST);
        const int abs_tol_status = abs_tol_.setData(abs_tol_data + offset, size_, memory::HOST);

        if (y_status != 0 || yp_status != 0 || f_status != 0 || abs_tol_status != 0)
        {
          Log::error() << "Component::bind - failed to bind vectors to system storage\n";
          return 1;
        }

        allocated_ = true;
        return 0;
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
      /**
       * @brief Allocate this component's state and residual vectors.
       */
      void allocateVectors(IdxT n)
      {

        y_.resize(n);
        yp_.resize(n);
        f_.resize(n);
        abs_tol_.resize(n);
      }

      /**
       * @brief Allocate this component's external variable vectors.
       */
      void allocateExternalVectors(IdxT n)
      {
        y_ext_.assign(static_cast<size_t>(n), ScalarT{});
        yp_ext_.assign(static_cast<size_t>(n), ScalarT{});
        variable_indices_ext_.assign(static_cast<size_t>(n), INVALID_INDEX<IdxT>);
      }

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

      IdxT              size_{0};
      IdxT              nnz_{0};
      /// Global (system-level) variable indices
      std::vector<IdxT> variable_indices_;
      std::vector<IdxT> variable_indices_ext_;
      /// Global (system-level) residual indices
      std::vector<IdxT> residual_indices_;
      std::vector<IdxT> residual_indices_ext_;

      VectorT              y_;
      std::vector<ScalarT> y_ext_;
      VectorT              yp_;
      std::vector<ScalarT> yp_ext_;
      std::vector<bool>    tag_;
      VectorT              abs_tol_;
      VectorT              f_;
      std::vector<ScalarT> f_ext_;
      bool                 allocated_{false};

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

      [[noreturn]] VectorT& yB() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const VectorT& yB() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] VectorT& ypB() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const VectorT& ypB() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] VectorT& param() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const VectorT& param() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] VectorT& param_up() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const VectorT& param_up() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] VectorT& param_lo() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const VectorT& param_lo() const override
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

      [[noreturn]] VectorT& getIntegrand() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const VectorT& getIntegrand() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] VectorT& getAdjointResidual() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const VectorT& getAdjointResidual() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] VectorT& getAdjointIntegrand() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const VectorT& getAdjointIntegrand() const override
      {
        throw NotImplementedError(__func__);
      }
    };

  } // namespace PhasorDynamics
} // namespace GridKit
