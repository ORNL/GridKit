#pragma once

#include <cassert>
#include <functional>
#include <memory>
#include <optional>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include <GridKit/Constants.hpp>
#include <GridKit/Model/PowerElectronics/Bus/MicrogridBus.hpp>
#include <GridKit/Model/PowerElectronics/CircuitComponent.hpp>
#include <GridKit/Model/PowerElectronics/CircuitNode.hpp>
#include <GridKit/Model/PowerElectronics/PartitionInterface/PartitionInterface.hpp>
#include <GridKit/Model/PowerElectronics/SystemModelPowerElectronics.hpp>
#include <GridKit/ScalarTraits.hpp>

namespace GridKit
{

  /**
   * @brief Represents a subset of a PowerElectronicsModel that can be evaluated
   *        independently.
   *
   * A SubsystemModel contains a collection of existing GridKit components and
   * nodes taken from a larger system. Variables owned by those components and
   * nodes become internal variables of the subsystem. Variables needed by those
   * components but owned outside the subsystem become external coupling
   * variables.
   *
   * Components normally store connection indices in the global system indexing.
   * During subsystem allocation, these indices are temporarily replaced with a
   * contiguous local subsystem indexing so that the subsystem can be evaluated
   * like an independent PowerElectronicsModel.
   *
   * External coupling values must be supplied before residual or Jacobian
   * evaluation, either directly through the external-data vectors or through a
   * forcing function.
   *
   * @tparam ScalarT Scalar type used by the model.
   * @tparam IdxT Index type used for variable and connection indices.
   */
  template <class ScalarT, typename IdxT>
  class SubsystemModel : public PowerElectronicsModel<ScalarT, IdxT>
  {
  public:
    struct ForcingData
    {
      std::vector<ScalarT> y;
      std::vector<ScalarT> yp;
    };

    using TimeFunction = std::function<ForcingData(ScalarT)>;

  protected:
    using SystemModel    = PowerElectronicsModel<ScalarT, IdxT>;
    using RealT          = typename CircuitComponent<ScalarT, IdxT>::RealT;
    using CsrMatrixT     = typename CircuitComponent<ScalarT, IdxT>::CsrMatrixT;
    using component_type = CircuitComponent<ScalarT, IdxT>;
    using node_type      = PowerElectronics::NodeBase<ScalarT, IdxT>;
    using interface_type = PartitionInterface<ScalarT, IdxT>;

    using SystemModel::abs_tol_;
    using SystemModel::allocated_;
    using SystemModel::allocateVectors;
    using SystemModel::alpha_;
    using SystemModel::connection_nodes_;
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

    SubsystemModel()
      : SystemModel(false)
    {
    }

    /**
     * @brief Default constructor for the system model
     *
     * @post System model parameters set as default
     */

    SubsystemModel(bool use_jac)
      : SystemModel(use_jac)
    {
    }

    ~SubsystemModel() override
    {
      for (auto* comp : interfaces_)
      {
        delete comp;
      }
      // SubsystemModel does not own components_/nodes_ — they belong to (and are
      // deleted by) the parent system model.
      components_.clear();
      nodes_.clear();
    }

    /**
     * @brief Allocate the subsystem for independent evaluation.
     *
     * Converts the selected components and nodes from global system indexing to
     * local subsystem indexing, allocates storage for the subsystem's internal
     * variables, and creates storage for coupling variables that lie outside the
     * subsystem.
     *
     * Component pointers are then redirected to either the subsystem's internal
     * vectors or its external coupling-data vectors. Finally, the subsystem
     * Jacobian structure is assembled using only entries whose row and column
     * belong to internal subsystem variables.
     *
     * @pre Components and nodes added to the subsystem must already be allocated
     *      by the parent system.
     *
     * @post Component and node connection indices use local subsystem indexing,
     *       and the subsystem is ready for residual and Jacobian evaluation.
     *
     * @return 0 on success, otherwise an error code returned during allocation.
     */
    int allocate() override
    {
      if (int err_code = buildIndexMappings())
      {
        return err_code;
      }

      if (int err_code = mapGlobalToLocal())
      {
        return err_code;
      }

      n_intern_ = internal_map_.size();
      n_extern_ = external_map_.size();
      size_     = n_intern_;

      // Allocate subsystem vectors.
      y_ext_data_.resize(n_extern_);
      yp_ext_data_.resize(n_extern_);
      f_ext_data_.resize(n_extern_);

      connection_nodes_ = std::make_unique<IdxT[]>(size_);
      if (!allocated_)
      {
        allocateVectors(static_cast<IdxT>(size_), true);
        abs_tol_.setToZero(memory::HOST);
      }

      tag_.resize(size_);

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
          }
        }
      }

      {
        // The offset for each component's internal variables in the system vector.
        // They start at 0, and are stacked on top of each other.
        size_t component_internal_idx = 0;
        for (component_type* comp : components_)
        {
          // Update component internal pointers to their correct offsets
          comp->setInternalPointer(&y_int_[component_internal_idx]);
          comp->setInternalDerivativePointer(&yp_int_[component_internal_idx]);
          comp->setInternalResidualPointer(&f_int_[component_internal_idx]);

          component_internal_idx += comp->getInternalSize();

          const auto& external_indices = comp->getExternIndices();

          for (IdxT local_index : external_indices)
          {
            const IdxT connection_index = comp->getNodeConnection(local_index);

            // A variable can be external from the component's point of view while still
            // being owned by this subsystem. In that case, connect the component directly
            // to the corresponding entry in the subsystem's internal vectors.
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

            // Otherwise the variable is owned outside this subsystem. Connect the
            // component to the subsystem's external coupling-data storage instead.
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

      // Allocation always rebuilds the system Jacobian and its COO-to-CSR map.
      delete csr_jac_;
      csr_jac_ = nullptr;

      delete[] map_to_csr_;
      map_to_csr_ = nullptr;

      // Evaluate component Jacobians to get sparsity
      for (component_type* component : components_)
      {
        component->evaluateJacobian();
      }

      // Check whether a Jacobian entry belongs to the subsystem Jacobian.
      // Only entries whose row and column are both internal subsystem
      // variables are retained.
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
     * @brief Update the subsystem external state and derivative data.
     *
     * If a forcing function is provided, evaluate it at the current subsystem
     * time and copy the returned coupling values into the external state vectors.
     *
     * @post y_ext_data_ and yp_ext_data_ contain the external coupling values
     *       returned by the forcing function, if one is set.
     *
     * @throws std::runtime_error If the forcing function returns vectors whose
     *         sizes do not match the subsystem external-data vectors.
     *
     * @return 0 on success.
     */
    int distributeExternalVectors()
    {

      if (forcing_function_)
      {
        const auto forcing = (*forcing_function_)(time_);

        if (forcing.y.size() != y_ext_data_.size() || forcing.yp.size() != yp_ext_data_.size())
        {
          throw std::runtime_error(
              "SubsystemModel::distributeExternalVectors: forcing function "
              "returned vectors with incorrect sizes.");
        }

        std::copy(forcing.y.begin(), forcing.y.end(), y_ext_data_.begin());
        std::copy(forcing.yp.begin(), forcing.yp.end(), yp_ext_data_.begin());
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
      if (int err_code = distributeExternalVectors())
      {
        return err_code;
      }

      return SystemModel::evaluateInternalResidual();
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
      if (int err_code = distributeExternalVectors())
      {
        return err_code;
      }

      return SystemModel::evaluateJacobian();
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
    void addComponent(component_type* component)
    {
      if (!component->isAllocated())
      {
        throw std::logic_error(
            "SubsystemModel::addComponent: cannot add an unallocated component.");
      }

      if (connections_are_local_)
      {
        throw std::logic_error(
            "SubsystemModel::addComponent: cannot add component while "
            "in local-indexed state. Call release() first.");
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
      if (!node->isAllocated())
      {
        throw std::logic_error(
            "SubsystemModel::addNode: cannot add an unallocated node.");
      }

      if (connections_are_local_)
      {
        throw std::logic_error(
            "SubsystemModel::addNode: cannot add node while in "
            "local-indexed state. Call release() first.");
      }

      SystemModel::addNode(node);
    }

    /**
     * @brief Add a partition interface to the subsystem.
     *
     * Adds the partition interface to the subsystem's component list and keeps a
     * separate reference to it in the interface list.
     *
     * @param component Pointer to the interface component to add.
     */
    void addInterface(interface_type* component)
    {
      addComponent(component);
      interfaces_.push_back(component);
    }

    /**
     * @brief Restore the subsystem topology to global system indexing.
     *
     * Reverses local subsystem connection indices with the original
     * global system indices. The local internal/external mappings are
     * then cleared so that components or nodes can safely be added or removed.
     *
     * @pre Component and node connections may use local subsystem indices.
     *
     * @post Component and node connections use their original global indices,
     *       the subsystem mappings are cleared, and the subsystem is marked
     *       unallocated.
     *
     * @return int 0 if successful, positive if there's a recoverable error, negative if unrecoverable
     */
    int release()
    {
      if (int err_code = mapLocalToGlobal())
      {
        return err_code;
      }

      internal_map_.clear();
      external_map_.clear();

      allocated_ = false;

      return 0;
    }

    const std::vector<IdxT>& getExternalDataIndices() const
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

    void setForcingFunction(TimeFunction function)
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
    /**
     * @brief Restore local subsystem connection indices to global system indices.
     *
     * Reverses \ref mapGlobalToLocal(). Internal subsystem indices are translated
     * through the subsystem connection-node table, while external subsystem
     * indices are translated through external_data_indices_.
     *
     * The component/node connectivity is unchanged; only the index representation
     * is restored.
     *
     * @pre Component and node connections use local subsystem indices.
     *
     * @post Component and node connections use their original global indices.
     *
     * @return 0 on success.
     */
    int mapLocalToGlobal()
    {
      if (!connections_are_local_)
      {
        return 0;
      }

      for (component_type* component : components_)
      {

        for (IdxT i = 0; i < component->size(); i++)
        {
          const IdxT index = component->getNodeConnection(i);

          if (index == neg1_)
          {
            continue;
          }
          else if (index < this->getInternalSize())
          {
            component->setConnectionNodes(i, this->getNodeConnection(index));
          }
          else
          {
            component->setConnectionNodes(i, external_data_indices_[index - this->getInternalSize()]);
          }
        }
      }

      for (node_type* node : nodes_)
      {

        for (IdxT i = 0; i < node->size(); i++)
        {
          const IdxT index = node->getNodeConnection(i).idx_;

          if (index == neg1_)
          {
            continue;
          }
          node->setConnectionNodes(i, this->getNodeConnection(index));
        }
      }

      connections_are_local_ = false;

      return 0;
    }

    /**
     * @brief Replace global connection indices with subsystem-local indices.
     *
     * Uses the mappings created by \ref buildIndexMappings() to rewrite the connection
     * indices stored by every component and node. Variables owned by the subsystem
     * use `internal_map_`; variables owned outside the subsystem use `external_map_`.
     *
     * This changes only the indexing used to identify connections; it does not
     * change the physical component/node connectivity.
     *
     * @pre `internal_map_` and `external_map_` have been constructed from global
     *      connection indices.
     *
     * @post All valid component and node connections use local subsystem indices.
     *
     * @return 0 on success.
     */
    int mapGlobalToLocal()
    {
      if (connections_are_local_)
      {
        return 0;
      }

      for (component_type* component : components_)
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

      for (node_type* node : nodes_)
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

      connections_are_local_ = true;

      return 0;
    }

    /**
     * @brief Build the global-to-local variable mappings for the subsystem.
     *
     * Examines the connection indices of all components and nodes in the
     * subsystem and divides the referenced variables into two groups:
     *
     * - Internal variables are owned by a component or node in this subsystem.
     * - External variables are required by a subsystem component but are owned
     *   outside the subsystem.
     *
     * Internal variables receive local indices first. External coupling variables
     * are then assigned indices immediately after the internal range. This gives
     * every variable referenced by the subsystem a unique local index while
     * preserving its original global index in the corresponding map.
     *
     * @pre Component and node connections use global system indices.
     *
     * @post internal_map_ and external_map_ contain the global-to-local mappings
     *       needed to convert the subsystem topology to local indexing.
     *
     * @return 0 on success.
     */
    int buildIndexMappings()
    {

      if (connections_are_local_)
      {
        return 0;
      }

      internal_map_.clear();
      external_map_.clear();

      size_t component_internal_idx = 0;
      // Pass 1: Add variables owned internally by subsystem components.
      for (component_type* comp : components_)
      {
        const auto& extern_indices = comp->getExternIndices();

        for (IdxT i = 0; i < comp->size(); i++)
        {
          const IdxT index = comp->getNodeConnection(i);

          if (index != neg1_ && !extern_indices.contains(i))
          {
            internal_map_[index] = component_internal_idx++;
          }
        }
      }

      // Pass 2: Add variables owned by subsystem nodes.
      for (node_type* node : nodes_)
      {

        for (IdxT i = 0; i < node->size(); i++)
        {
          const IdxT index = node->getNodeConnection(i).idx_;

          if (index != neg1_)
          {
            internal_map_[index] = component_internal_idx++;
          }
        }
      }

      // Pass 3: Add component dependencies that are owned outside the subsystem.
      size_t component_external_idx = component_internal_idx;
      for (component_type* comp : components_)
      {
        auto extern_indices = comp->getExternIndices();
        for (IdxT j = 0; j < comp->size(); j++)
        {
          if (!extern_indices.contains(j))
          {
            continue;
          }

          const IdxT index = comp->getNodeConnection(j);

          if (index != neg1_ && !internal_map_.contains(index) && !external_map_.contains(index))
          {
            external_map_[index] = component_external_idx++;
          }
        }
      }

      return 0;
    }

    /**
     *@brief Maps global system indices to local subsystem indices for internal variables.
     */
    std::unordered_map<IdxT, IdxT> internal_map_;

    /**
     * @brief Maps global system indices to local subsystem indices for external variables.
     */
    std::unordered_map<IdxT, IdxT> external_map_;

    /**
     * @brief Global system index corresponding to each entry in the external subsystem vectors.
     */
    std::vector<IdxT> external_data_indices_;

    /**
     * @brief subsystem external State, derivative, and residual vectors.
     */
    std::vector<ScalarT> y_ext_data_;

    /**
     * @brief subsystem external derivative
     */
    std::vector<ScalarT> yp_ext_data_;

    /**
     * @brief subsystem external derivative
     */
    std::vector<ScalarT> f_ext_data_;

    /**
     * @brief Optional forcing function used to provide external subsystem data.
     *
     * It is continuous function that can be sampled at different times
     */
    std::optional<TimeFunction> forcing_function_;

    /**
     * @brief Partition interfaces owned by the subsystem.
     *
     * These interfaces expose the contributions of components participating
     * in the partition split.
     */
    std::vector<interface_type*> interfaces_;

    /**
     * @brief Keeps track of whether the components are in local state or global state
     */
    bool connections_are_local_{false};

  }; // class SubsystemModel

} // namespace GridKit