#pragma once

#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/CommonMath.hpp>
#include <GridKit/Constants.hpp>
#include <GridKit/Model/EMT/SignalNode/SignalNode.hpp>
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
     * submodels. A component owns its internal variables and one residual row
     * per internal variable, reads external variables and their derivatives
     * through bound signal nodes, and accumulates external residual
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
      using SignalT    = SignalNode<ScalarT, IdxT>;
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
       * The component's own rows occupy the first ownSize() entries of the
       * slice; each registered submodel is bound to the sub-slice that
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

        const IdxT own_size = ownSize();

        const int y_status       = y_.setData(y_data + offset, own_size, memory::HOST);
        const int yp_status      = yp_.setData(yp_data + offset, own_size, memory::HOST);
        const int f_status       = f_.setData(f_data + offset, own_size, memory::HOST);
        const int abs_tol_status = abs_tol_.setData(abs_tol_data + offset, own_size, memory::HOST);

        if (y_status != 0 || yp_status != 0 || f_status != 0 || abs_tol_status != 0)
        {
          Log::error() << "Component::bind - failed to bind vectors to system storage\n";
          return 1;
        }

        for (size_t i = 0; i < submodels_.size(); ++i)
        {
          const int submodel_status = submodels_[i]->bind(y, yp, f, abs_tol, offset + submodel_offsets_[i]);
          if (submodel_status != 0)
          {
            return submodel_status;
          }
        }

        allocated_ = true;
        return 0;
      }

      /**
       * @brief Assign the global index of one local variable.
       *
       * Indices for rows owned by a registered submodel are forwarded into
       * the submodel's local index space in addition to being recorded here.
       */
      int setVariableIndex(IdxT local_index, IdxT global_index)
      {
        variable_indices_[static_cast<size_t>(local_index)] = global_index;
        routeSubmodelIndex(local_index, global_index, true);
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
       * Indices for rows owned by a registered submodel are forwarded into
       * the submodel's local index space in addition to being recorded here.
       */
      int setResidualIndex(IdxT local_index, IdxT global_index)
      {
        residual_indices_[static_cast<size_t>(local_index)] = global_index;
        routeSubmodelIndex(local_index, global_index, false);
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
       * @brief Register the signal node backing one external variable slot.
       */
      int setExternalVariableNode(IdxT slot, SignalT* node)
      {
        external_variable_nodes_[static_cast<size_t>(slot)] = node;
        return 0;
      }

      /**
       * @brief Register the signal node whose residual row receives one
       * external residual contribution.
       */
      int setExternalResidualNode(IdxT row, SignalT* node)
      {
        external_residual_nodes_[static_cast<size_t>(row)] = node;
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
        for (auto* submodel : submodels_)
        {
          submodel->updateTime(t, a);
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
       * next evaluateJacobian() call. Registered submodels are invalidated
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

        for (auto* submodel : submodels_)
        {
          submodel->resetJacobianStructure();
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
       * Equal to size() for a component without submodels.
       */
      IdxT ownSize() const
      {
        if (submodels_.empty())
        {
          return size_;
        }
        return own_size_;
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
        external_variable_nodes_.assign(static_cast<size_t>(n_vars), nullptr);
        f_ext_.assign(static_cast<size_t>(n_rows), ScalarT{});
        residual_indices_ext_.assign(static_cast<size_t>(n_rows), INVALID_INDEX<IdxT>);
        external_residual_nodes_.assign(static_cast<size_t>(n_rows), nullptr);
      }

      /**
       * @brief Bind a three-phase port over three consecutive local variables.
       *
       * Each phase node exposes the variable, its derivative, and its
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
          port.nodes[static_cast<size_t>(n)].set(&y[local_first + n],
                                                 &yp[local_first + n],
                                                 &f[local_first + n],
                                                 &getVariableIndex(local_first + n),
                                                 &getResidualIndex(local_first + n));
        }
        return 0;
      }

      /**
       * @brief Gather external variables, derivatives, and index maps through
       * the registered signal nodes.
       *
       * Models with latched defaults for optional signals override this
       * method, call the base implementation, and patch the unattached slots.
       */
      virtual void gatherExternalVariables()
      {
        for (size_t i = 0; i < external_variable_nodes_.size(); ++i)
        {
          auto* node = external_variable_nodes_[i];
          if (node != nullptr && node->linked())
          {
            y_ext_[i]                = node->read();
            variable_indices_ext_[i] = node->getVariableIndex();
            yp_ext_[i]               = ScalarT{};
            if (node->derivativeLinked())
            {
              yp_ext_[i] = node->readDerivative();
            }
          }
        }

        for (size_t i = 0; i < external_residual_nodes_.size(); ++i)
        {
          auto* node = external_residual_nodes_[i];
          if (node != nullptr && node->residualLinked())
          {
            residual_indices_ext_[i] = node->getResidualIndex();
          }
        }
      }

      /**
       * @brief Accumulate this component's external residual contributions
       * into the residual rows owned by the connected components.
       */
      void scatterExternalResidual()
      {
        for (size_t i = 0; i < external_residual_nodes_.size(); ++i)
        {
          auto* node = external_residual_nodes_[i];
          if (node != nullptr && node->residualLinked())
          {
            node->accumulateResidual(f_ext_[i]);
          }
        }
      }

      /**
       * @brief Register a submodel owning the next rows of this component.
       *
       * The caller adds the submodel's size into size_ before allocation.
       */
      int registerSubmodel(Component* submodel)
      {
        submodel_offsets_.push_back(own_size_ + submodelSize());
        submodels_.push_back(submodel);
        return 0;
      }

      int allocateSubmodels()
      {
        for (auto* submodel : submodels_)
        {
          const int status = submodel->allocate();
          if (status != 0)
          {
            return status;
          }
        }
        return 0;
      }

      int initializeSubmodels()
      {
        for (auto* submodel : submodels_)
        {
          const int status = submodel->initialize();
          if (status != 0)
          {
            return status;
          }
        }
        return 0;
      }

      int tagDifferentiableSubmodels()
      {
        for (size_t i = 0; i < submodels_.size(); ++i)
        {
          submodels_[i]->tagDifferentiable();
          const auto& submodel_tag = submodels_[i]->tag();
          const auto  offset       = static_cast<size_t>(submodel_offsets_[i]);
          for (size_t j = 0; j < submodel_tag.size(); ++j)
          {
            tag_[offset + j] = submodel_tag[j];
          }
        }
        return 0;
      }

      int setAbsoluteToleranceSubmodels(RealT rel_tol)
      {
        for (auto* submodel : submodels_)
        {
          submodel->setAbsoluteTolerance(rel_tol);
        }
        return 0;
      }

      int evaluateSubmodelInternalResiduals()
      {
        for (auto* submodel : submodels_)
        {
          const int status = submodel->evaluateInternalResidual();
          if (status != 0)
          {
            return status;
          }
        }
        return 0;
      }

      int evaluateSubmodelExternalResiduals()
      {
        for (auto* submodel : submodels_)
        {
          const int status = submodel->evaluateExternalResidual();
          if (status != 0)
          {
            return status;
          }
        }
        return 0;
      }

      int evaluateSubmodelJacobians()
      {
        for (auto* submodel : submodels_)
        {
          const int status = submodel->evaluateJacobian();
          if (status != 0)
          {
            return status;
          }
        }
        return 0;
      }

      /**
       * @brief Copy submodel Jacobian triplets into this component's buffers.
       *
       * Submodel triplets already carry global indices because index
       * assignment was routed through setVariableIndex and setResidualIndex,
       * so this component presents one COO to the system.
       *
       * @pre This component's Jacobian buffers have capacity for its own
       * entries plus every submodel's entries.
       */
      int appendSubmodelJacobians()
      {
        for (auto* submodel : submodels_)
        {
          auto* coo = submodel->getCooJacobian();
          if (coo == nullptr)
          {
            continue;
          }
          const IdxT  submodel_nnz = coo->getNnz();
          const auto* rows         = coo->getRowData();
          const auto* cols         = coo->getColData();
          const auto* vals         = coo->getValues();
          for (IdxT i = 0; i < submodel_nnz; ++i)
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
      /// Rows owned directly by this component, excluding submodel rows
      IdxT              own_size_{0};
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

      /// Signal node backing each external variable slot
      std::vector<SignalT*> external_variable_nodes_;
      /// Signal node whose residual row receives each external residual row
      std::vector<SignalT*> external_residual_nodes_;

      /// Registered submodels, in row order after this component's own rows
      std::vector<Component*> submodels_;
      /// Local row offset of each registered submodel
      std::vector<IdxT>       submodel_offsets_;

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
      IdxT submodelSize() const
      {
        IdxT total = 0;
        for (const auto* submodel : submodels_)
        {
          total += const_cast<Component*>(submodel)->size();
        }
        return total;
      }

      void routeSubmodelIndex(IdxT local_index, IdxT global_index, bool is_variable)
      {
        if (submodels_.empty() || local_index < own_size_)
        {
          return;
        }
        for (size_t i = submodels_.size(); i > 0; --i)
        {
          const auto submodel_i = i - 1;
          if (local_index >= submodel_offsets_[submodel_i])
          {
            const IdxT submodel_local = local_index - submodel_offsets_[submodel_i];
            if (is_variable)
            {
              submodels_[submodel_i]->setVariableIndex(submodel_local, global_index);
            }
            else
            {
              submodels_[submodel_i]->setResidualIndex(submodel_local, global_index);
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
