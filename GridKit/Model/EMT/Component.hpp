#pragma once

#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/CommonMath.hpp>
#include <GridKit/Constants.hpp>
#include <GridKit/Model/EMT/Signal/Signal.hpp>
#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Utilities/Errors.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace EMT
  {
    using Log = ::GridKit::Utilities::Logger;

    /**
     * @brief Component model implementation base class.
     *
     * Every EMT model is a Component, including buses and models that own
     * operators. A component owns its internal variables and one residual row
     * per internal variable, reads external variables and their derivatives
     * through bound signals, and accumulates external residual
     * contributions into residual rows owned by other components.
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
      using SignalT    = Signal<ScalarT, IdxT>;
      using Port3T     = Port3<ScalarT, IdxT>;

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

      /**
       * @brief Stable ordering key for component initialization.
       *
       * Containers initialize all leaves in ascending order. Components with
       * the same value retain their hierarchy order.
       */
      virtual int initializationOrder() const noexcept
      {
        return 1;
      }

      virtual int evaluateInternalResidual()
      {
        return this->evaluateResidual();
      }

      virtual int evaluateExternalResidual()
      {
        return 0;
      }

      IdxT size() override
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
       * The component's own rows occupy the first equationSize() entries of the
       * slice; each embedded operator is bound to the sub-slice that
       * follows, in registration order.
       *
       * After binding, the component reads and writes system storage directly
       * and allocate() will not allocate local vector data. Rebinding is
       * allowed and refreshes the aliases, e.g. after system storage is
       * reallocated when the topology changes.
       *
       * Only HOST data is bound because EMT currently evaluates models on the
       * CPU. Supporting DEVICE data would also require sharing the matching
       * DEVICE pointer and keeping the HOST and DEVICE copies in sync. This
       * bind operation does neither, so DEVICE access is unsupported.
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
      virtual int bind(VectorT& y, VectorT& yp, VectorT& f, VectorT& abs_tol, IdxT offset)
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

        if (size_ == 0)
        {
          allocated_ = true;
          return 0;
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

        const IdxT own_size = equationSize();

        const int y_status       = y_.setData(y_data + offset, own_size, memory::HOST);
        const int yp_status      = yp_.setData(yp_data + offset, own_size, memory::HOST);
        const int f_status       = f_.setData(f_data + offset, own_size, memory::HOST);
        const int abs_tol_status = abs_tol_.setData(abs_tol_data + offset, own_size, memory::HOST);

        if (y_status != 0 || yp_status != 0 || f_status != 0 || abs_tol_status != 0)
        {
          Log::error() << "Component::bind - failed to bind vectors to system storage\n";
          return 1;
        }

        for (size_t i = 0; i < operators_.size(); ++i)
        {
          const int operator_status = operators_[i]->bind(y, yp, f, abs_tol, offset + operator_offsets_[i]);
          if (operator_status != 0)
          {
            return operator_status;
          }
        }

        allocated_ = true;
        return 0;
      }

      /**
       * @brief Assign the global index of one local variable.
       *
       * Indices for rows owned by an embedded operator are forwarded into
       * the operator's local index space in addition to being recorded here.
       */
      int setVariableIndex(IdxT local_index, IdxT global_index)
      {
        variable_indices_[static_cast<size_t>(local_index)] = global_index;
        routeOperatorIndex(local_index, global_index, true);
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

      /**
       * @brief Assign the global index of one local residual row.
       *
       * Indices for rows owned by an embedded operator are forwarded into
       * the operator's local index space in addition to being recorded here.
       */
      int setResidualIndex(IdxT local_index, IdxT global_index)
      {
        residual_indices_[static_cast<size_t>(local_index)] = global_index;
        routeOperatorIndex(local_index, global_index, false);
        return 0;
      }

      /**
       * @brief Assign contiguous system indices to this equation block.
       *
       * Containers override this once for an entire subtree. The default
       * implementation also preserves the private operator rows used by the
       * rational EMT devices.
       */
      virtual int assignGlobalIndices(IdxT first)
      {
        for (IdxT j = 0; j < size_; ++j)
        {
          setVariableIndex(j, first + j);
          setResidualIndex(j, first + j);
        }
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

      /**
       * @brief Register the signal backing one external variable slot.
       */
      int setExternalVariableSignal(IdxT slot, SignalT* signal)
      {
        external_variable_signals_[static_cast<size_t>(slot)] = signal;
        return 0;
      }

      const std::vector<SignalT*>& externalVariableSignals() const
      {
        return external_variable_signals_;
      }

      bool hasComputedInputs() const
      {
        for (const auto* signal : external_variable_signals_)
        {
          if (signal != nullptr && signal->computed())
          {
            return true;
          }
        }
        return false;
      }

      /// Maximum number of global columns represented by one input signal.
      size_t externalJacobianExpansion() const
      {
        size_t                      expansion = 1;
        typename SignalT::GradientT gradient;
        for (const auto* signal : external_variable_signals_)
        {
          if (signal != nullptr && signal->computed())
          {
            gradient.clear();
            signal->appendGradient(gradient);
            expansion = std::max(expansion, gradient.size());
          }
        }
        return expansion;
      }

      /**
       * @brief Register the signal whose residual row receives one
       * external residual contribution.
       */
      int setExternalResidualSignal(IdxT row, SignalT* signal)
      {
        external_residual_signals_[static_cast<size_t>(row)] = signal;
        return 0;
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
        for (auto* op : operators_)
        {
          op->updateTime(t, a);
        }
      }

      virtual int setGridKitComponentID(IdxT) = 0;

      IdxT getGridKitComponentID() const
      {
        return gridkit_component_id_;
      }

      /**
       * @brief Invalidate this component's Jacobian structure.
       *
       * Discrete state changes can change which Jacobian entries are exactly
       * zero, and exact zeros are dropped from the discovered sparsity
       * pattern, so the cached COO and CSR structures must be rebuilt by the
       * next evaluateJacobian() call. Registered operators are invalidated
       * along with their parent.
       */
      virtual void resetJacobianStructure()
      {
        delete coo_jac_;
        coo_jac_ = nullptr;

        delete csr_jac_;
        csr_jac_ = nullptr;

        delete[] map_to_csr_;
        map_to_csr_ = nullptr;

        nnz_ = 0;

        for (auto* op : operators_)
        {
          op->resetJacobianStructure();
        }
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
       * @brief Number of rows owned directly by this component.
       *
       * Equal to size() for a component without operators.
       */
      IdxT equationSize() const
      {
        if (operators_.empty())
        {
          return size_;
        }
        return equation_size_;
      }

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
       * @brief Allocate this component's external variable and residual vectors.
       *
       * @param[in] n_vars - Number of external variable slots
       * @param[in] n_rows - Number of external residual rows
       */
      void allocateExternalVectors(IdxT n_vars, IdxT n_rows)
      {
        y_ext_.assign(static_cast<size_t>(n_vars), ScalarT{});
        yp_ext_.assign(static_cast<size_t>(n_vars), ScalarT{});
        variable_indices_ext_.assign(static_cast<size_t>(n_vars), INVALID_INDEX<IdxT>);
        external_variable_signals_.assign(static_cast<size_t>(n_vars), nullptr);
        f_ext_.assign(static_cast<size_t>(n_rows), ScalarT{});
        residual_indices_ext_.assign(static_cast<size_t>(n_rows), INVALID_INDEX<IdxT>);
        external_residual_signals_.assign(static_cast<size_t>(n_rows), nullptr);
      }

      /**
       * @brief Bind a three-phase port over three consecutive local variables.
       *
       * Each phase signal exposes the variable, its derivative, and its
       * residual row along with the global indices, so the port can serve as
       * a connection surface for values, derivative reads, and residual
       * accumulation.
       *
       * @pre The component vectors are allocated or bound, and the index
       * vectors are sized.
       */
      int bindPort(Port3T& port, IdxT local_first)
      {
        auto* y  = y_.getData();
        auto* yp = yp_.getData();
        auto* f  = f_.getData();
        for (IdxT n = 0; n < 3; ++n)
        {
          port.signals[static_cast<size_t>(n)].set(&y[local_first + n],
                                                   &yp[local_first + n],
                                                   &f[local_first + n],
                                                   &getVariableIndex(local_first + n),
                                                   &getResidualIndex(local_first + n));
        }
        return 0;
      }

      /**
       * @brief Gather external variables, derivatives, and index maps through
       * the registered signals.
       *
       * Models with latched defaults for optional signals override this
       * method, call the base implementation, and patch the unattached slots.
       */
      virtual void gatherExternalVariables()
      {
        for (size_t i = 0; i < external_variable_signals_.size(); ++i)
        {
          auto* signal = external_variable_signals_[i];
          if (signal != nullptr && signal->linked())
          {
            y_ext_[i]                = signal->read();
            variable_indices_ext_[i] = signal->getVariableIndex();
            yp_ext_[i]               = ScalarT{};
            if (signal->derivativeLinked())
            {
              yp_ext_[i] = signal->readDerivative();
            }
          }
        }

        for (size_t i = 0; i < external_residual_signals_.size(); ++i)
        {
          auto* signal = external_residual_signals_[i];
          if (signal != nullptr && signal->residualLinked())
          {
            residual_indices_ext_[i] = signal->getResidualIndex();
          }
        }
      }

      /**
       * @brief Accumulate this component's external residual contributions
       * into the residual rows owned by the connected components.
       */
      void scatterExternalResidual()
      {
        for (size_t i = 0; i < external_residual_signals_.size(); ++i)
        {
          auto* signal = external_residual_signals_[i];
          if (signal != nullptr && signal->residualLinked())
          {
            signal->accumulateResidual(f_ext_[i]);
          }
        }
      }

      /**
       * @brief Add an embedded operator owning the next rows of this component.
       *
       * The caller adds the operator's size into size_ before allocation.
       *
       * @pre equation_size_ is already set, because the operator's local row
       * offset is recorded here.
       */
      int addOperator(Component* op)
      {
        operator_offsets_.push_back(equation_size_ + operatorSize());
        operators_.push_back(op);
        return 0;
      }

      int allocateOperators()
      {
        for (auto* op : operators_)
        {
          const int status = op->allocate();
          if (status != 0)
          {
            return status;
          }
        }
        return 0;
      }

      int initializeOperators()
      {
        for (auto* op : operators_)
        {
          const int status = op->initialize();
          if (status != 0)
          {
            return status;
          }
        }
        return 0;
      }

      int tagDifferentiableOperators()
      {
        for (size_t i = 0; i < operators_.size(); ++i)
        {
          operators_[i]->tagDifferentiable();
          const auto& operator_tag = operators_[i]->tag();
          const auto  offset       = static_cast<size_t>(operator_offsets_[i]);
          for (size_t j = 0; j < operator_tag.size(); ++j)
          {
            tag_[offset + j] = operator_tag[j];
          }
        }
        return 0;
      }

      int setAbsoluteToleranceOperators(RealT rel_tol)
      {
        for (auto* op : operators_)
        {
          op->setAbsoluteTolerance(rel_tol);
        }
        return 0;
      }

      int evaluateOperatorInternalResiduals()
      {
        for (auto* op : operators_)
        {
          const int status = op->evaluateInternalResidual();
          if (status != 0)
          {
            return status;
          }
        }
        return 0;
      }

      int evaluateOperatorExternalResiduals()
      {
        for (auto* op : operators_)
        {
          const int status = op->evaluateExternalResidual();
          if (status != 0)
          {
            return status;
          }
        }
        return 0;
      }

      int evaluateOperatorJacobians()
      {
        for (auto* op : operators_)
        {
          const int status = op->evaluateJacobian();
          if (status != 0)
          {
            return status;
          }
        }
        return 0;
      }

      /**
       * @brief Copy operator Jacobian triplets into this component's buffers.
       *
       * Operator triplets already carry global indices because index
       * assignment was routed through setVariableIndex and setResidualIndex,
       * so this component presents one COO to the system.
       *
       * @pre This component's Jacobian buffers have capacity for its own
       * entries plus every operator's entries.
       */
      int appendOperatorJacobians()
      {
        for (auto* op : operators_)
        {
          auto* coo = op->getCooJacobian();
          if (coo == nullptr)
          {
            continue;
          }
          const IdxT  operator_nnz = coo->getNnz();
          const auto* rows         = coo->getRowData();
          const auto* cols         = coo->getColData();
          const auto* vals         = coo->getValues();
          for (IdxT i = 0; i < operator_nnz; ++i)
          {
            J_rows_buffer_[nnz_] = rows[i];
            J_cols_buffer_[nnz_] = cols[i];
            J_vals_buffer_[nnz_] = vals[i];
            ++nnz_;
          }
        }
        return 0;
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
      /// Rows owned directly by this component, excluding operator rows
      IdxT              equation_size_{0};
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

      /// Signal backing each external variable slot
      std::vector<SignalT*> external_variable_signals_;
      /// Signal whose residual row receives each external residual row
      std::vector<SignalT*> external_residual_signals_;

      /// Registered operators, in row order after this component's own rows
      std::vector<Component*> operators_;
      /// Local row offset of each embedded operator
      std::vector<IdxT>       operator_offsets_;

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

      // A nonzero default keeps derivative entries in the Jacobian pattern
      // when sparsity is discovered before the first solver step.
      RealT time_{0.0};
      RealT alpha_{1.0};

      using NotImplementedError = GridKit::Utilities::NotImplementedError;

    private:
      IdxT operatorSize() const
      {
        IdxT total = 0;
        for (const auto* op : operators_)
        {
          total += const_cast<Component*>(op)->size();
        }
        return total;
      }

      void routeOperatorIndex(IdxT local_index, IdxT global_index, bool is_variable)
      {
        if (operators_.empty() || local_index < equation_size_)
        {
          return;
        }
        for (size_t i = operators_.size(); i > 0; --i)
        {
          const auto operator_i = i - 1;
          if (local_index >= operator_offsets_[operator_i])
          {
            const IdxT operator_local = local_index - operator_offsets_[operator_i];
            if (is_variable)
            {
              operators_[operator_i]->setVariableIndex(operator_local, global_index);
            }
            else
            {
              operators_[operator_i]->setResidualIndex(operator_local, global_index);
            }
            return;
          }
        }
      }

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

  } // namespace EMT
} // namespace GridKit
