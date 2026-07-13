

#pragma once

#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

#include <GridKit/Constants.hpp>
#include <GridKit/Model/PowerElectronics/CircuitComponent.hpp>
#include <GridKit/Model/PowerElectronics/CircuitNode.hpp>
#include <GridKit/Model/PowerElectronics/NodeBase.hpp>
#include <GridKit/ScalarTraits.hpp>

namespace GridKit
{
  template <class ScalarT, typename IdxT>
  class PowerElectronicsModel : public CircuitComponent<ScalarT, IdxT>
  {
    using RealT          = typename CircuitComponent<ScalarT, IdxT>::RealT;
    using CsrMatrixT     = typename CircuitComponent<ScalarT, IdxT>::CsrMatrixT;
    using component_type = CircuitComponent<ScalarT, IdxT>;
    using node_type      = PowerElectronics::NodeBase<ScalarT, IdxT>;

    using CircuitComponent<ScalarT, IdxT>::size_;
    using CircuitComponent<ScalarT, IdxT>::n_intern_;
    using CircuitComponent<ScalarT, IdxT>::n_extern_;
    using CircuitComponent<ScalarT, IdxT>::nnz_;
    using CircuitComponent<ScalarT, IdxT>::time_;
    using CircuitComponent<ScalarT, IdxT>::alpha_;
    using CircuitComponent<ScalarT, IdxT>::y_;
    using CircuitComponent<ScalarT, IdxT>::y_int_;
    using CircuitComponent<ScalarT, IdxT>::yp_;
    using CircuitComponent<ScalarT, IdxT>::yp_int_;
    using CircuitComponent<ScalarT, IdxT>::f_;
    using CircuitComponent<ScalarT, IdxT>::f_int_;
    using CircuitComponent<ScalarT, IdxT>::tag_;
    using CircuitComponent<ScalarT, IdxT>::abs_tol_;
    using CircuitComponent<ScalarT, IdxT>::allocated_;
    using CircuitComponent<ScalarT, IdxT>::allocateVectors;

  public:
    /**
     * @brief Default constructor for the system model
     *
     * @post System model parameters set as default
     */
    PowerElectronicsModel()
    {
      // By default don't use the jacobian
      use_jac_ = false;
    }

    /**
     * @brief Constructor for the system model
     *
     * @param[in] use_jac Boolean to choose if to use jacobian
     *
     * @post System model parameters set as input
     */
    PowerElectronicsModel(bool use_jac = false)
    {
      // Can choose if to use jacobian
      use_jac_ = use_jac;
    }

    /**
     * @brief Destructor for the system model
     *
     * @pre System components are allocated
     *
     * @post System components are deallocated
     *
     */
    virtual ~PowerElectronicsModel()
    {
      for (auto comp : this->components_)
      {
        delete comp;
      }
      delete csr_jac_;
      delete[] map_to_csr_;
    }

    /**
     * @brief Will check if each component has jacobian avalible. If one doesn't have it, return false.
     *
     *
     * @return true if all components have jacobian
     * @return false otherwise
     */
    bool hasJacobian() final
    {
      if (!this->use_jac_)
        return false;

      for (const auto& component : components_)
      {
        if (!component->hasJacobian())
        {
          return false;
        }
      }
      return true;
    }

    /**
     * @brief Allocate system vectors and construct the system CSR Jacobian
     *
     * @param[in] s size of the vector (total number of unknowns)
     *
     * @post System model vectors allocated with size s
     * @post CSR Jacobian sparsity pattern is computed
     * @post COO->CSR mapping is computed
     * @post Every component's \ref CircuitComponent::y_int_, \ref CircuitComponent::yp_int_, and \ref CircuitComponent::f_int_ pointers
     * are set to their appropriate offsets in the system vector, allowing them to directly access their internal variables, derivatives,
     * and residuals.
     *
     * @return int 0 if successful, positive if there's a recoverable error, negative if unrecoverable
     */
    int allocate() final
    {
      size_t component_internal_size = 0;
      for (component_type* comp : components_)
      {
        component_internal_size += comp->getInternalSize();
      }

      size_t node_internal_size = 0;
      for (node_type* node : nodes_)
      {
        node_internal_size += node->getInternalSize();
      }

      n_intern_ = component_internal_size + node_internal_size;
      n_extern_ = 0;
      size_     = n_intern_ + n_extern_;

      if (!allocated_)
      {
        allocateVectors(static_cast<IdxT>(size_));
      }

      tag_.resize(size_);

      { // Start node internal indexing after all component internals for proper KLU ordering
        size_t node_internal_idx = component_internal_size;
        for (node_type* node : nodes_)
        {
          node->allocate();

          for (size_t i = 0; i < node->getInternalSize(); i++)
          {
            node->setExternalConnectionNodes(i, node_internal_idx);
            node_internal_idx++;
          }
        }
      }

      {
        // The offset for each component's internal variables in the system vector.
        // They start at 0, and are stacked on top of each other.
        size_t component_internal_idx = 0;
        auto*  y                      = y_.getData();
        auto*  yp                     = yp_.getData();
        auto*  f                      = f_.getData();
        for (component_type* comp : components_)
        {
          comp->allocate();

          // Update component internal pointers to their correct offsets
          comp->setInternalPointer(&y[component_internal_idx]);
          comp->setInternalDerivativePointer(&yp[component_internal_idx]);
          comp->setInternalResidualPointer(&f[component_internal_idx]);

          const auto& external_indices = comp->getExternIndices();
          for (size_t i = 0; i < comp->size(); i++)
          {
            if (!external_indices.contains(i))
            {
              comp->setExternalConnectionNodes(i, component_internal_idx);
              component_internal_idx++;
            }
          }
        }
      }

      // Evaluate component Jacobians to get sparsity
      distributeVectors();
      for (component_type* component : components_)
      {
        component->evaluateJacobian();
      }

      // Count the number of non-zeros
      IdxT nnz_dup = 0;
      for (const component_type* component : components_)
      {
        const IdxT* r   = component->jacobianCooRows();
        const IdxT* c   = component->jacobianCooCols();
        IdxT        nnz = component->nnz();

        for (IdxT i = 0; i < nnz; ++i)
        {
          if (component->getNodeConnection(r[i]) != neg1_ && component->getNodeConnection(c[i]) != neg1_)
          {
            ++nnz_dup;
          }
        }
      }

      // Allocate COO triplet arrays (we own these until we hand off to CsrMatrix)
      IdxT*  rows_dup = new IdxT[nnz_dup];
      IdxT*  cols_dup = new IdxT[nnz_dup];
      RealT* vals_dup = new RealT[nnz_dup];

      IdxT counter = 0;
      for (const auto& component : components_)
      {
        const IdxT*  r   = component->jacobianCooRows();
        const IdxT*  c   = component->jacobianCooCols();
        const RealT* v   = component->jacobianCooValues();
        IdxT         nnz = component->nnz();

        for (IdxT i = 0; i < nnz; ++i)
        {
          if (component->getNodeConnection(r[i]) != neg1_ && component->getNodeConnection(c[i]) != neg1_)
          {
            rows_dup[counter] = component->getNodeConnection(r[i]);
            cols_dup[counter] = component->getNodeConnection(c[i]);
            vals_dup[counter] = v[i];
            counter++;
          }
        }
      }

      // Build the system COO Jacobian
      LinearAlgebra::CooMatrix<RealT, IdxT> jac(size_, size_, nnz_dup, &rows_dup, &cols_dup, &vals_dup);

      // Populate CSR data with sort and deduplicate
      IdxT* row_ptrs = jac.getCsrRowData();

      // Deduplicated nnz
      nnz_ = jac.getNnz();

      // Allocate cols/vals with deduplicated nnz
      IdxT*  cols = new IdxT[nnz_];
      RealT* vals = new RealT[nnz_];

      std::copy(jac.getColData(), jac.getColData() + nnz_, cols);
      std::copy(jac.getValues(), jac.getValues() + nnz_, vals);

      // Create the CSR Jacobian
      csr_jac_ = new CsrMatrixT(size_, size_, nnz_, &row_ptrs, &cols, &vals);

      const IdxT* map_to_sorted = jac.getMapToSorted();
      const IdxT* map_to_dedup  = jac.getMapToDeduplicated();

      // Build a mappping from original COO index to CSR index
      map_to_csr_ = new IdxT[nnz_dup];
      for (IdxT i = 0; i < nnz_dup; ++i)
      {
        map_to_csr_[map_to_sorted[i]] = map_to_dedup[i];
      }

      return 0;
    }

    /**
     * @brief Set intial y and y' of each component
     *
     * @return int 0 if successful, positive if there's a recoverable error, negative if unrecoverable
     */
    int initialize() final
    {
      // Initialize components
      for (const auto& component : components_)
      {
        component->initialize();
      }
      this->distributeVectors();

      return 0;
    }

    /**
     * @brief Distribute y and y' to each component based of node connection graph
     *
     * @post Each component has y and y' set
     *
     * @return int 0 if successful, positive if there's a recoverable error, negative if unrecoverable
     */
    int distributeVectors()
    {
      auto* y_system  = y_.getData();
      auto* yp_system = yp_.getData();

      for (component_type* component : components_)
      {
        auto*                   y         = component->y().getData();
        auto*                   yp        = component->yp().getData();
        const std::set<size_t>& externals = component->getExternIndices();

        for (size_t j : externals)
        {
          if (component->getNodeConnection(j) != neg1_)
          {
            y[j]  = y_system[component->getNodeConnection(j)];
            yp[j] = yp_system[component->getNodeConnection(j)];
          }
          else
          {
            y[j]  = 0.0;
            yp[j] = 0.0;
          }
        }
      }
      return 0;
    }

    int tagDifferentiable() final
    {
      return 0;
    }

    /**
     * @brief Compute the absolute tolerance for each variable in the model
     *
     * @param rel_tol The relative tolerance which can be used to pick the
     *        absolute tolerance.
     * @tparam ScalarT Scalar data type
     * @tparam IdxT Index data type
     * @return int 0 if successful, non-zero otherwise.
     *
     * This represents a "noise" level close to zero for which pure relative
     * error cannot be used.
     */
    int setAbsoluteTolerance(RealT rel_tol) final
    {
      std::fill(abs_tol_.getData(), abs_tol_.getData() + abs_tol_.getSize(), rel_tol);
      return 0;
    }

    /**
     * @brief Evaluate Residuals at each component then collect them
     *
     * @return int 0 if successful, positive if there's a recoverable error, negative if unrecoverable
     */
    int evaluateInternalResidual() final
    {
      auto* f = f_.getData();

      for (IdxT i = 0; i < this->f_.getSize(); i++)
      {
        f[i] = 0.0;
      }

      this->distributeVectors();

      // Update system residual vector

      // Evaluate component internal residuals - this is embarassingly parallel
      for (component_type* component : components_)
      {
        if (int err_code = component->evaluateInternalResidual())
          return err_code;
      }

      for (component_type* component : components_)
      {
        if (int err_code = component->evaluateExternalResidual())
          return err_code;

        auto*                   residual  = component->getResidual().getData();
        const std::set<size_t>& externals = component->getExternIndices();

        for (size_t j : externals)
        {
          //@todo should do a different grounding check
          if (component->getNodeConnection(j) != neg1_)
          {
            f[component->getNodeConnection(j)] += residual[j];
          }
        }
      }

      return 0;
    }

    /**
     * @todo implement this for nested systems
     */
    int evaluateExternalResidual() final
    {
      return 0;
    }

    /**
     * @brief Creates the system Jacobian representing \alpha dF/dy' + dF/dy
     *
     * Updates the CSR Jacobian values using the per-component mappings
     * computed during allocate().
     *
     * @return int 0 if successful, positive if there's a recoverable error, negative if unrecoverable
     */
    int evaluateJacobian() final
    {
      distributeVectors();

      // Zero out values
      RealT* vals = csr_jac_->getValues();
      for (IdxT i = 0; i < csr_jac_->getNnz(); ++i)
      {
        vals[i] = 0.0;
      }

      // Update CSR values from component Jacobians
      IdxT counter = 0;
      for (const auto& component : components_)
      {
        component->evaluateJacobian();

        const IdxT*  r   = component->jacobianCooRows();
        const IdxT*  c   = component->jacobianCooCols();
        const RealT* v   = component->jacobianCooValues();
        IdxT         nnz = component->nnz();

        for (IdxT i = 0; i < nnz; ++i)
        {
          if (component->getNodeConnection(r[i]) != neg1_ && component->getNodeConnection(c[i]) != neg1_)
          {
            vals[map_to_csr_[counter]] += v[i];
            ++counter;
          }
        }
      }

      jac_call_count_++;
      return 0;
    }

    /**
     * @brief Evaluate integrands for the system quadratures.
     */
    int evaluateIntegrand() final
    {
      return 0;
    }

    /**
     * @brief Initialize system adjoint.
     *
     * Updates variables and optimization parameters, then initializes
     * adjoints locally and copies them to the system adjoint vector.
     */
    int initializeAdjoint() final
    {
      return 0;
    }

    /**
     * @brief Compute adjoint residual for the system model.
     *
     *
     */
    int evaluateAdjointResidual() final
    {
      return 0;
    }

    /**
     * @brief Evaluate adjoint integrand for the system model.
     *
     *
     */
    int evaluateAdjointIntegrand() final
    {
      return 0;
    }

    /**
     * @brief Distribute time and time scaling for each component
     *
     * @param t
     * @param a
     */
    void updateTime(RealT t, RealT a) final
    {
      for (const auto& component : components_)
      {
        component->updateTime(t, a);
      }
      time_  = t;
      alpha_ = a;
    }

    CsrMatrixT* getCsrJacobian() const override
    {
      return csr_jac_;
    }

    void addComponent(component_type* component)
    {
      components_.push_back(component);
    }

    void addNode(node_type* node)
    {
      nodes_.push_back(node);
    }

  private:
    static constexpr IdxT neg1_ = INVALID_INDEX<IdxT>;

    std::vector<component_type*> components_;
    std::vector<node_type*>      nodes_;

    IdxT*       map_to_csr_{nullptr};
    CsrMatrixT* csr_jac_{nullptr};

    int  jac_call_count_{0};
    bool use_jac_;

  }; // class PowerElectronicsModel

} // namespace GridKit
