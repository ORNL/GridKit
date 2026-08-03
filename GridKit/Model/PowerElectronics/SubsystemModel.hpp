

#pragma once

#include <cassert>
#include <functional>
#include <optional>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <GridKit/Constants.hpp>
#include <GridKit/Model/PowerElectronics/Bus/MicrogridBus.hpp>
#include <GridKit/Model/PowerElectronics/CircuitComponent.hpp>
#include <GridKit/Model/PowerElectronics/CircuitNode.hpp>
#include <GridKit/Model/PowerElectronics/SystemModelPowerElectronics.hpp>
#include <GridKit/ScalarTraits.hpp>

namespace GridKit
{

  template <class ScalarT, typename IdxT>
  class SubsystemModel : public PowerElectronicsModel<ScalarT, IdxT>
  {
  public:
    using ForcingData  = std::tuple<std::vector<ScalarT>, std::vector<ScalarT>>;
    using TimeFunction = std::function<ForcingData(ScalarT)>;

  protected:
    using SystemModel    = PowerElectronicsModel<ScalarT, IdxT>;
    using RealT          = typename CircuitComponent<ScalarT, IdxT>::RealT;
    using CsrMatrixT     = typename CircuitComponent<ScalarT, IdxT>::CsrMatrixT;
    using component_type = CircuitComponent<ScalarT, IdxT>;
    using node_type      = PowerElectronics::NodeBase<ScalarT, IdxT>;

    using SystemModel::abs_tol_;
    using SystemModel::allocated_;
    using SystemModel::allocateVectors;
    using SystemModel::alpha_;
    using SystemModel::f_int_;
    using SystemModel::n_extern_;
    using SystemModel::n_intern_;
    using SystemModel::nnz_;
    using SystemModel::size_;
    using SystemModel::tag_;
    using SystemModel::time_;
    using SystemModel::y_int_;
    using SystemModel::yp_int_;

    using SystemModel::components_;
    using SystemModel::csr_jac_;
    using SystemModel::jac_call_count_;
    using SystemModel::map_to_csr_;
    using SystemModel::neg1_;
    using SystemModel::nodes_;
    using SystemModel::use_jac_;

  public:
    /**
     * @brief Default constructor for the system model
     *
     * @post System model parameters set as default
     */

    SubsystemModel(bool use_jac = false)
      : SystemModel(use_jac)
    {
    }

    ~SubsystemModel() override
    {
      // SubsystemModel does not own components_/nodes_ — they belong to (and are
      // deleted by) the parent system model.
      components_.clear();
      nodes_.clear();
    }

    int allocate() override
    {

      hold();

      n_intern_ = internal_map_.size();
      n_extern_ = external_map_.size();
      size_     = n_intern_;

      if (!allocated_)
      {
        allocateVectors(static_cast<IdxT>(size_), true);
      }

      CircuitComponent<ScalarT, IdxT>::allocate();

      // Allocation always rebuilds the system Jacobian and its COO-to-CSR map.
      delete csr_jac_;
      csr_jac_ = nullptr;

      delete[] map_to_csr_;
      map_to_csr_ = nullptr;

      tag_.resize(size_);

      // Allocate subsystem vectors.
      y_ext_data_.resize(n_extern_);
      yp_ext_data_.resize(n_extern_);
      f_ext_data_.resize(n_extern_);

      external_data_indices_.resize(n_extern_);

      // Store the mapping from local subsystem indices back to their global system indices
      for (const auto [global_idx, local_idx] : internal_map_)
      {
        this->setConnectionNodes(local_idx, global_idx);
      }

      // Store the global indices of all external coupling variables.
      for (const auto [global_idx, local_idx] : external_map_)
      {
        external_data_indices_[local_idx - n_intern_] = global_idx;
      }

      tag_.resize(size_);

      { // Start node internal indexing after all component internals for proper KLU ordering
        size_t node_internal_idx;
        for (node_type* node : nodes_)
        {
          for (size_t i = 0; i < node->getInternalSize(); i++)
          {
            node_internal_idx = node->getNodeConnection(i).idx_;

            ExternalConnection<ScalarT, IdxT> node_connection{
                .y_   = y_int_ + node_internal_idx,
                .yp_  = yp_int_ + node_internal_idx,
                .f_   = f_int_ + node_internal_idx,
                .idx_ = static_cast<IdxT>(node_internal_idx)};

            node->setExternalConnectionNodes(i, node_connection);
            node_internal_idx++;
          }
        }
      }

      {
        // The offset for each component's internal variables in the system vector.
        // They start at 0, and are stacked on top of each other.
        size_t component_internal_idx = 0;
        for (component_type* comp : components_)
        {
          if (!comp->isAllocated())
          {
            comp->allocate();
          }
          // Update component internal pointers to their correct offsets
          comp->setInternalPointer(&y_int_[component_internal_idx]);
          comp->setInternalDerivativePointer(&yp_int_[component_internal_idx]);
          comp->setInternalResidualPointer(&f_int_[component_internal_idx]);

          component_internal_idx += comp->getInternalSize();

          const auto& external_indices = comp->getExternIndices();

          for (IdxT local_index : external_indices)
          {
            const IdxT connection_index = comp->getNodeConnection(local_index);

            // External to the component but internal to the subsystem.
            if (connection_index < n_intern_)
            {
              ExternalConnection<ScalarT, IdxT> connection{
                  .y_   = y_int_ + connection_index,
                  .yp_  = yp_int_ + connection_index,
                  .f_   = f_int_ + connection_index,
                  .idx_ = connection_index};

              comp->setExternalConnectionNodes(local_index, connection);

              continue;
            }

            // External to both the component and subsystem.
            const IdxT external_offset = connection_index - static_cast<IdxT>(n_intern_);

            ExternalConnection<ScalarT, IdxT> connection{
                .y_   = &y_ext_data_[external_offset],
                .yp_  = &yp_ext_data_[external_offset],
                .f_   = &f_ext_data_[external_offset],
                .idx_ = connection_index};

            comp->setExternalConnectionNodes(local_index, connection);
          }
        }
      }

      // Evaluate component Jacobians to get sparsity
      for (component_type* component : components_)
      {
        component->evaluateJacobian();
      }

      auto isValidEntry = [this](IdxT row, IdxT col)
      {
        if (row == neg1_ || col == neg1_)
        {
          return false;
        }

        const bool row_is_internal = row < this->getInternalSize();
        const bool col_is_internal = col < this->getInternalSize();

        return (row_is_internal && col_is_internal);
      };

      IdxT nnz_dup = 0;

      for (const component_type* component : components_)
      {
        const IdxT* r   = component->jacobianCooRows();
        const IdxT* c   = component->jacobianCooCols();
        const IdxT  nnz = component->nnz();

        for (IdxT i = 0; i < nnz; ++i)
        {
          const IdxT row = component->getNodeConnection(r[i]);
          const IdxT col = component->getNodeConnection(c[i]);

          if (isValidEntry(row, col))
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

      for (const component_type* component : components_)
      {
        const IdxT*  r   = component->jacobianCooRows();
        const IdxT*  c   = component->jacobianCooCols();
        const RealT* v   = component->jacobianCooValues();
        const IdxT   nnz = component->nnz();

        for (IdxT i = 0; i < nnz; ++i)
        {
          const IdxT row = component->getNodeConnection(r[i]);
          const IdxT col = component->getNodeConnection(c[i]);

          if (!isValidEntry(row, col))
          {
            continue;
          }

          rows_dup[counter] = row;
          cols_dup[counter] = col;
          vals_dup[counter] = v[i];

          ++counter;
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

      allocated_ = true;
      return 0;
    }

    /**
     * @brief Distribute y and y' to each component based of node connection graph
     *
     * @post Each component has y and y' set
     *
     * @return int 0 if successful, positive if there's a recoverable error, negative if unrecoverable
     */
    int distributeExternalVectors()
    {

      if (forcing_function_)
      {
        const auto [y_forcing, yp_forcing] = (*forcing_function_)(time_);

        assert(y_forcing.size() == y_ext_data_.size());
        assert(yp_forcing.size() == yp_ext_data_.size());

        std::copy(y_forcing.begin(), y_forcing.end(), y_ext_data_.begin());
        std::copy(yp_forcing.begin(), yp_forcing.end(), yp_ext_data_.begin());
      }

      return 0;
    }

    /**
     * @brief Evaluate Residuals at each component then collect them
     *
     * @return int 0 if successful, positive if there's a recoverable error, negative if unrecoverable
     */
    int evaluateInternalResidual() override
    {
      distributeExternalVectors();

      SystemModel::evaluateInternalResidual();

      return 0;
    }

    /**
     * @brief Creates the system Jacobian representing \f$\alpha dF/dy' + dF/dy\f$
     *
     * Updates the CSR Jacobian values using the per-component mappings
     * computed during allocate().
     *
     * @return int 0 if successful, positive if there's a recoverable error, negative if unrecoverable
     */
    int evaluateJacobian() override
    {
      distributeExternalVectors();

      SystemModel::evaluateJacobian();

      return 0;
    }

    /**
     * @brief Add a component to the subsystem.
     *
     * Rejected while connections are in the local-indexed state, since the
     * component's stored connection indices would otherwise be interpreted
     * inconsistently with the rest of the subsystem. Call release() first.
     *
     * @param[in] component Component to add.
     */
    void
    addComponent(component_type* component)

    {
      if (is_comp_in_local_state)
      {
        std::cerr << "SubsystemModel::addComponent: cannot add component while "
                     "in local-indexed state. Call release() first.\n";
        return;
      }

      SystemModel::addComponent(component);
    }

    /**
     * @brief Add a node to the subsystem.
     *
     * Rejected while connections are in the local-indexed state, since the
     * node's stored connection indices would otherwise be interpreted
     * inconsistently with the rest of the subsystem. Call release() first.
     *
     * @param[in] node Node to add.
     */
    void addNode(node_type* node)
    {
      if (is_comp_in_local_state)
      {
        std::cerr << "SubsystemModel::addNode: cannot add node while in "
                     "local-indexed state. Call release() first.\n";
        return;
      }

      SystemModel::addNode(node);
    }

    /**
     * @brief Acquire the subsystem's global-to-local index mappings.
     *
     * Builds the internal/external index maps via buildIndexMappings(),
     * then converts component/node connections from global to local
     * indices via mapGlobalToLocal().
     *
     * Call this after modifying subsystem topology (adding or removing
     * components/nodes) and before allocate(), to (re)establish local
     * indexing. Pairs with release(), which undoes this.
     *
     * @return int 0 if successful, positive if there's a recoverable error, negative if unrecoverable
     */
    int hold()
    {
      buildIndexMappings();

      if (int err_code = mapGlobalToLocal())
      {
        return err_code;
      }

      return 0;
    }

    // /**
    //  * @brief Release the subsystem's global-to-local index mappings.
    //  *
    //  * If components/nodes currently hold local indices, they are first
    //  * converted back to their original global indices via
    //  * mapLocalToGlobal(). The internal/external index maps built by
    //  * buildIndexMappings() are then cleared, and the model is marked as
    //  * not allocated.
    //  *
    //  *
    //  * Call this before modifying subsystem topology (adding or removing
    //  * components/nodes), so stale local indices aren't left dangling.
    //  *
    //  * @return int 0 if successful, positive if there's a recoverable error, negative if unrecoverable
    //  */
    // int release()
    // {
    //   if (int err_code = mapLocalToGlobal())
    //   {
    //     return err_code;
    //   }

    //   internal_map_.clear();
    //   external_map_.clear();

    //   allocated_ = false;

    //   return 0;
    // }

    const std::vector<IdxT>& getExternalDataIndices()
    {
      return external_data_indices_;
    }

    std::vector<ScalarT>& getExternalDataY()
    {
      return y_ext_data_;
    }

    std::vector<ScalarT>& getExternalDataYP()
    {
      return yp_ext_data_;
    }

    std::vector<ScalarT>& getExternalDataF()
    {
      return f_ext_data_;
    }

    void setTimeFunction(TimeFunction function)
    {
      forcing_function_ = std::move(function);
    }

    const std::unordered_map<IdxT, IdxT>& getInternalMap() const
    {
      return internal_map_;
    }

    const std::unordered_map<IdxT, IdxT>& getExternalMap() const
    {
      return external_map_;
    }

  private:
    // /**
    //  * @brief Replace local subsystem connection indices with their original global indices
    //  * for every component and node.
    //  *
    //  * Inverse of mapGlobalToLocal(). Internal connections (index < internal
    //  * size) are mapped back through this subsystem's own node-connection
    //  * table; external connections are mapped back through
    //  * @c external_data_indices_. Entries equal to @c neg1_ represent no
    //  * connection and are left unchanged.
    //  *
    //  * @return Always returns 0.
    //  */
    // int mapLocalToGlobal()
    // {
    //   if (!is_comp_in_local_state)
    //   {
    //     return 0;
    //   }

    //   for (auto* component : components_)
    //   {

    //     for (IdxT i = 0; i < component->size(); i++)
    //     {
    //       IdxT index = component->getNodeConnection(i);

    //       if (index == neg1_)
    //       {
    //         continue;
    //       }
    //       else if (index < this->getInternalSize())
    //       {
    //         component->setExternalConnectionNodes(i, this->getNodeConnection(index));
    //       }
    //       else
    //       {
    //         component->setExternalConnectionNodes(i, external_data_indices_[index - this->getInternalSize()]);
    //       }
    //     }
    //   }

    //   for (auto* node : nodes_)
    //   {

    //     for (IdxT i = 0; i < node->size(); i++)
    //     {
    //       IdxT index = node->getNodeConnection(i);

    //       if (index == neg1_)
    //       {
    //         continue;
    //       }
    //       node->setExternalConnectionNodes(i, this->getNodeConnection(index));
    //     }
    //   }

    //   is_comp_in_local_state = false;

    //   return 0;
    // }

    /**
     * @brief Replace global connection indices with local subsystem indices for every component
     * and node.
     *
     * After the internal and external maps have been constructed, each component
     * and node still stores the original global connection indices. This function
     * traverses every connection and replaces each global index with its
     * corresponding local subsystem index.
     *
     * Internal connections are mapped using @c internal_map_, while connections
     * that reference variables outside the subsystem are mapped using
     * @c external_map_. Entries equal to @c neg1_ represent no connection and
     * are left unchanged.
     *
     * @return Always returns 0.
     */

    int mapGlobalToLocal()
    {
      if (is_comp_in_local_state)
      {
        return 0;
      }

      for (auto* component : components_)
      {
        for (IdxT i = 0; i < component->size(); ++i)
        {
          const IdxT index = component->getNodeConnection(i);

          if (index == neg1_)
          {
            continue;
          }

          if (internal_map_.contains(index))
          {
            component->setConnectionNodes(i, internal_map_.at(index));
          }
          else
          {
            component->setConnectionNodes(i, external_map_.at(index));
          }
        }
      }

      for (auto* node : nodes_)
      {
        for (IdxT i = 0; i < node->size(); ++i)
        {
          const IdxT index = node->getNodeConnection(i).idx_;

          if (index == neg1_)
          {
            continue;
          }

          node->setConnectionNodes(i, internal_map_.at(index));
        }
      }

      is_comp_in_local_state = true;

      return 0;
    }

    /**
     * @brief Construct the mappings from global variable indices to local
     * subsystem indices.
     *
     * The subsystem stores its variables in a contiguous local ordering. This
     * function assigns each global variable index a corresponding local index and
     * separates them into internal and external variables.
     *
     * Internal variables are those owned by this subsystem, while external
     * variables are owned by neighboring subsystems but are required for residual
     * and Jacobian evaluation.
     */
    void buildIndexMappings()
    {

      if (is_comp_in_local_state)
      {
        return;
      }

      internal_map_.clear();
      external_map_.clear();

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
          IdxT index = node->getNodeConnection(i).idx_;

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

    std::unordered_map<IdxT, IdxT> internal_map_;
    std::unordered_map<IdxT, IdxT> external_map_;
    std::vector<IdxT>              external_data_indices_;

    std::vector<ScalarT> y_ext_data_;
    std::vector<ScalarT> yp_ext_data_;
    std::vector<ScalarT> f_ext_data_;

    std::optional<TimeFunction> forcing_function_;

    bool is_comp_in_local_state{false};

  }; // class SubsystemModel

} // namespace GridKit