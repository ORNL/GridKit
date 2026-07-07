

#pragma once

#include <cassert>
#include <string>
#include <unordered_map>
#include <vector>

#include <GridKit/Constants.hpp>
#include <GridKit/Model/PowerElectronics/Bus/MicrogridBus.hpp>
#include <GridKit/Model/PowerElectronics/CircuitComponent.hpp>
#include <GridKit/Model/PowerElectronics/CircuitNode.hpp>
#include <GridKit/ScalarTraits.hpp>

namespace GridKit
{
  template <class ScalarT, typename IdxT>
  class SubsystemModel : public CircuitComponent<ScalarT, IdxT>
  {
    using RealT          = typename CircuitComponent<ScalarT, IdxT>::RealT;
    using MatrixT        = typename CircuitComponent<ScalarT, IdxT>::MatrixT;
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
    using CircuitComponent<ScalarT, IdxT>::rel_tol_;
    using CircuitComponent<ScalarT, IdxT>::abs_tol_;
    using CircuitComponent<ScalarT, IdxT>::extern_indices_;

  public:
    /**
     * @brief Default constructor for the system model
     *
     * @post System model parameters set as default
     */

    SubsystemModel(bool root = true)
    {
      // Set system model parameters as default
      rel_tol_         = 1e-4;
      abs_tol_         = 1e-4;
      this->max_steps_ = 2000;
      // By default don't use the jacobian
      use_jac_         = false;
      root_            = root;
    }

    virtual ~SubsystemModel()
    {
    }

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

    int allocate() override
    {

      this->createGlobalToInternalMap();

      this->mapGlobalToLocal();

      n_intern_ = internal_map_.size();
      n_extern_ = external_map_.size();
      size_     = n_intern_;

      CircuitComponent<ScalarT, IdxT>::allocate();

      y_.resize(n_intern_);
      yp_.resize(n_intern_);
      f_.resize(n_intern_);

      // Allocate subsystem vectors.
      y_ext_.resize(n_extern_);
      yp_ext_.resize(n_extern_);
      f_ext_.resize(n_extern_);
      external_indices_.resize(n_extern_);

      // Store the mapping from local subsystem indices back to their global system indices
      for (const auto [global_idx, local_idx] : internal_map_)
      {
        this->setExternalConnectionNodes(local_idx, global_idx);
      }

      // Store the global indices of all external coupling variables.
      for (const auto [global_idx, local_idx] : external_map_)
      {
        external_indices_[local_idx - n_intern_] = global_idx;
      }

      
        size_t component_internal_idx = 0;
        for (component_type* comp : components_)
        {

          comp->setInternalPointer(&y_[component_internal_idx]);
          comp->setInternalDerivativePointer(&yp_[component_internal_idx]);
          comp->setInternalResidualPointer(&f_[component_internal_idx]);

          component_internal_idx += comp->getInternalSize();
        }
      

      // // Evaluate component Jacobians to get sparsity
      // distributeVectors();
      // for (component_type* component : components_)
      // {
      //   component->evaluateJacobian();
      // }

      // // Count the number of non-zeros
      // IdxT nnz_dup = 0;
      // for (const component_type* component : components_)
      // {
      //   const IdxT* r   = component->jacobianCooRows();
      //   const IdxT* c   = component->jacobianCooCols();
      //   IdxT        nnz = component->nnz();

      //   for (IdxT i = 0; i < nnz; ++i)
      //   {
      //     if(component->getNodeConnection(r[i])>= this->getInternalSize() && component->getNodeConnection(c[i])>= this->getInternalSize())
      //     {
      //       continue;
      //     }
      //     if (component->getNodeConnection(r[i]) != neg1_ && component->getNodeConnection(c[i]) != neg1_)
      //     {
      //       ++nnz_dup;
      //     }
      //   }
      // }

      // // Allocate COO triplet arrays (we own these until we hand off to CsrMatrix)
      // IdxT*  rows_dup = new IdxT[nnz_dup];
      // IdxT*  cols_dup = new IdxT[nnz_dup];
      // RealT* vals_dup = new RealT[nnz_dup];

      // IdxT counter = 0;
      // for (const auto& component : components_)
      // {
      //   const IdxT*  r   = component->jacobianCooRows();
      //   const IdxT*  c   = component->jacobianCooCols();
      //   const RealT* v   = component->jacobianCooValues();
      //   IdxT         nnz = component->nnz();

      //   for (IdxT i = 0; i < nnz; ++i)
      //   {
      //     if(component->getNodeConnection(r[i])>= this->getInternalSize() && component->getNodeConnection(c[i])>= this->getInternalSize())
      //     {
      //       continue;
      //     }

      //     if (component->getNodeConnection(r[i]) != neg1_ && component->getNodeConnection(c[i]) != neg1_)
      //     {
      //       rows_dup[counter] = component->getNodeConnection(r[i]);
      //       cols_dup[counter] = component->getNodeConnection(c[i]);
      //       vals_dup[counter] = v[i];
      //       counter++;
      //     }
      //   }
      // }

      // // Build the system COO Jacobian
      // LinearAlgebra::CooMatrix<RealT, IdxT> jac(size_, size_, nnz_dup, &rows_dup, &cols_dup, &vals_dup);

      // // Populate CSR data with sort and deduplicate
      // IdxT* row_ptrs = jac.getCsrRowData();

      // // Deduplicated nnz
      // nnz_ = jac.getNnz();

      // // Allocate cols/vals with deduplicated nnz
      // IdxT*  cols = new IdxT[nnz_];
      // RealT* vals = new RealT[nnz_];

      // std::copy(jac.getColData(), jac.getColData() + nnz_, cols);
      // std::copy(jac.getValues(), jac.getValues() + nnz_, vals);

      // // Create the CSR Jacobian
      // csr_jac_ = new CsrMatrixT(size_, size_, nnz_, &row_ptrs, &cols, &vals);

      // const IdxT* map_to_sorted = jac.getMapToSorted();
      // const IdxT* map_to_dedup  = jac.getMapToDeduplicated();

      // // Build a mappping from original COO index to CSR index
      // map_to_csr_ = new IdxT[nnz_dup];
      // for (IdxT i = 0; i < nnz_dup; ++i)
      // {
      //   map_to_csr_[map_to_sorted[i]] = map_to_dedup[i];
      // }

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
      for (component_type* component : components_)
      {
        std::vector<ScalarT>&   y         = component->y();
        std::vector<ScalarT>&   yp        = component->yp();
        const std::set<size_t>& externals = component->getExternIndices();

        for (size_t j : externals)
        {
          if (component->getNodeConnection(j) == neg1_)
          {
            y[j]  = 0.0;
            yp[j] = 0.0;
          }
          else if (component->getNodeConnection(j) < this->getInternalSize())
          {
            y[j]  = y_[component->getNodeConnection(j)];
            yp[j] = yp_[component->getNodeConnection(j)];
          }
          else
          {
            y[j]  = y_ext_[component->getNodeConnection(j) - this->getInternalSize()];
            yp[j] = yp_ext_[component->getNodeConnection(j) - this->getInternalSize()];
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
     * @brief Evaluate Residuals at each component then collect them
     *
     * @return int 0 if successful, positive if there's a recoverable error, negative if unrecoverable
     */
    int evaluateInternalResidual() final
    {
      for (IdxT i = 0; i < this->getInternalSize(); i++)
      {
        f_[i] = 0.0;
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

        const std::vector<ScalarT>& residual  = component->getResidual();
        const std::set<size_t>&     externals = component->getExternIndices();

        for (size_t j : externals)
        {
          //@todo should do a different grounding check
          if (component->getNodeConnection(j) != neg1_ && component->getNodeConnection(j) < this->getInternalSize())
          {
            f_[component->getNodeConnection(j)] += residual[j];
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

          if(component->getNodeConnection(r[i]) >= this->getInternalSize() && component->getNodeConnection(c[i]) >= this->getInternalSize())
          {
            continue;
          }

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

    /**
     * @brief print the system residual in COO format
     *
     * @param[in] filename
     * @param[in] title
     */
    void printResidualMatrixMarket(std::string filename, std::string title)
    {
      writeVectorToMatrixMarket(f_, filename, title);
    }

    /**
     * @brief print the system Jacobian in COO format
     *
     * @param[in] filename
     * @param[in] title
     */
    void printJacobianMatrixMarket(std::string filename, std::string title)
    {
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

    int mapGlobalToLocal()
    {
      for (auto* component : components_)
      {

        for (IdxT i = 0; i < component->size(); i++)
        {
          IdxT index = component->getNodeConnection(i);

          if (index == neg1_)
          {
            continue;
          }
          else if (internal_map_.count(index) > 0)
          {
            component->setExternalConnectionNodes(i, internal_map_.at(index));
          }
          else
          {
            component->setExternalConnectionNodes(i, external_map_.at(index));
          }
        }
      }

      for (auto* node : nodes_)
      {

        for (IdxT i = 0; i < node->size(); i++)
        {
          IdxT index = node->getNodeConnection(i);

          if (index == neg1_)
          {
            continue;
          }
          else if (internal_map_.count(index) > 0)
          {
            node->setExternalConnectionNodes(i, internal_map_.at(index));
          }
          else
          {
            node->setExternalConnectionNodes(i, external_map_.at(index));
          }
        }
      }

      return 0;
    }

    void createGlobalToInternalMap()
    {

      size_t component_internal_idx = 0;
      // First pass: Map global component indices to local subsystem indices.
      for (component_type* comp : components_)
      {
        for (IdxT i = 0; i < comp->size(); i++)
        {
          IdxT index = comp->getNodeConnection(i);

          auto extern_indices = comp->getExternIndices();

          if (index != neg1_ && !extern_indices.contains(i))
          {
            internal_map_[index] = component_internal_idx++;
          }
        }
      }

      for (node_type* node : nodes_)
      {

        for (IdxT i = 0; i < node->size(); i++)
        {
          IdxT index = node->getNodeConnection(i);

          if (index != neg1_)
          {
            internal_map_[index] = component_internal_idx++;
          }
        }
      }

      for (component_type* comp : components_)
      {
        auto extern_indices = comp->getExternIndices();
        for (IdxT j = 0; j < comp->size(); j++)
        {
          if (!extern_indices.contains(j))
          {
            continue;
          }
          IdxT index = comp->getNodeConnection(j);

          if (internal_map_.count(index) < 1 && external_map_.count(index) < 1 && index != neg1_)
          {
            external_map_[index] = component_internal_idx++;
          }
        }
      }
    }

    std::vector<IdxT>& getExternalIndices()
    {
      return external_indices_;
    }

    std::vector<ScalarT>& getExternalDataY()
    {
      return y_ext_;
    }

    std::vector<ScalarT>& getExternalDataYP()
    {
      return yp_ext_;
    }

    std::vector<ScalarT>& getExternalDataF()
    {
      return f_ext_;
    }

  private:
    static constexpr IdxT neg1_ = INVALID_INDEX<IdxT>;

    std::vector<component_type*>   components_;
    std::vector<node_type*>        nodes_;
    std::unordered_map<IdxT, IdxT> internal_map_;
    std::unordered_map<IdxT, IdxT> external_map_;
    std::vector<IdxT>              external_indices_;
    bool                           root_;
    std::vector<ScalarT>           y_ext_;
    std::vector<ScalarT>           yp_ext_;
    std::vector<ScalarT>           f_ext_;

    IdxT*       map_to_csr_{nullptr};
    CsrMatrixT* csr_jac_{nullptr};

    int  jac_call_count_{0};
    bool use_jac_;

  }; // class SubsystemModel

} // namespace GridKit
