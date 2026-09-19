/**
 * @file BusPartitionInterface.hpp
 * @author Abdourahman Barry (abdourahman@vt.edu)
 *
 */
#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
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
   * @todo Find a better name for this class and its base class.
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
    using VectorT        = typename SystemModel::VectorT;
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
      if (this->isAllocated())
      {
        restoreGlobalConnection();
      }

      if (int err_code = buildConnectionMap())
      {
        return err_code;
      }

      // Allocate subsystem vectors.
      y_ext_data_.resize(n_extern_);
      yp_ext_data_.resize(n_extern_);
      f_ext_data_.resize(n_extern_);

      connection_nodes_ = std::make_unique<IdxT[]>(n_intern_ + n_extern_);
      if (!this->isAllocated())
      {
        allocateVectors(static_cast<IdxT>(size_), true);
        abs_tol_.setToZero(memory::HOST);
      }

      tag_.resize(size_);

      external_data_indices_.resize(n_extern_);

      // Store the mapping from local subsystem indices back to their global system indices
      for (const auto [global_idx, local_idx] : global_to_local_map_)
      {
        this->setConnectionNodes(local_idx, global_idx);
      }

      // Store the global indices of all external coupling variables.
      size_t counter = 0;
      for (size_t i = n_intern_; i < (n_intern_ + n_extern_); i++)
      {
        external_data_indices_[counter] = this->getNodeConnection(i);
        counter++;
      }

      // Start node internal indexing after all component internals for proper KLU ordering
      for (node_type* node : nodes_)
      {
        for (size_t i = 0; i < node->getInternalSize(); i++)
        {
          const IdxT node_global_connection = node->getNodeConnection(i).idx_;
          const IdxT node_local_connection  = global_to_local_map_.at(node_global_connection);

          ExternalConnection<ScalarT, IdxT> node_connection{
              .y_   = y_int_ + node_local_connection,
              .yp_  = yp_int_ + node_local_connection,
              .f_   = f_int_ + node_local_connection,
              .idx_ = static_cast<IdxT>(node_local_connection)};

          node->setExternalConnectionNodes(i, node_connection);
        }
      }

      ScalarT* y_ext_data_ptr  = y_ext_data_.getData();
      ScalarT* yp_ext_data_ptr = yp_ext_data_.getData();
      ScalarT* f_ext_data_ptr  = f_ext_data_.getData();

      size_t component_internal_idx = 0;
      for (component_type* comp : components_)
      {
        // Update component internal pointers to their correct offsets
        comp->setInternalPointer(&y_int_[component_internal_idx]);
        comp->setInternalDerivativePointer(&yp_int_[component_internal_idx]);
        comp->setInternalResidualPointer(&f_int_[component_internal_idx]);

        component_internal_idx += comp->getInternalSize();

        const auto& external_indices = comp->getExternIndices();

        for (IdxT i = 0; i < comp->size(); i++)
        {
          const IdxT comp_global_connection = comp->getNodeConnection(i);
          const IdxT comp_local_connection  = global_to_local_map_.at(comp_global_connection);

          // Internal component variables use their subsystem-local connection indices.
          if (!external_indices.contains(i))
          {
            comp->setConnectionNodes(i, comp_local_connection);
            continue;
          }

          if (comp_local_connection < this->getInternalSize())
          {
            // This variable is external to the component but internal to the subsystem,
            // since it is owned by another component within the subsystem.

            ExternalConnection<ScalarT, IdxT> connection{
                .y_   = y_int_ + comp_local_connection,
                .yp_  = yp_int_ + comp_local_connection,
                .f_   = f_int_ + comp_local_connection,
                .idx_ = comp_local_connection};

            comp->setExternalConnectionNodes(i, connection);
          }
          else
          {
            // Otherwise the variable is owned outside this subsystem. Connect the
            // component to the subsystem's external coupling-data storage instead.
            const IdxT external_offset = comp_local_connection - static_cast<IdxT>(this->getInternalSize());

            ExternalConnection<ScalarT, IdxT> connection{
                .y_   = y_ext_data_ptr + external_offset,
                .yp_  = yp_ext_data_ptr + external_offset,
                .f_   = f_ext_data_ptr + external_offset,
                .idx_ = comp_local_connection};

            comp->setExternalConnectionNodes(i, connection);
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
        if (row == INVALID_INDEX<IdxT> || col == INVALID_INDEX<IdxT>)
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

      counter = 0;
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

        if (forcing.y.size() != y_ext_data_.getSize() || forcing.yp.size() != yp_ext_data_.getSize())
        {
          throw std::runtime_error(
              "SubsystemModel::distributeExternalVectors: forcing function "
              "returned vectors with incorrect sizes.");
        }

        std::copy(forcing.y.begin(), forcing.y.end(), y_ext_data_.getData());
        std::copy(forcing.yp.begin(), forcing.yp.end(), yp_ext_data_.getData());
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
     * @brief Restores the subsystem to global system indexing.
     *
     * Restores the original global connection indices of all components and nodes,
     * then clears the subsystem's global-to-local connection map. The subsystem is
     * marked unallocated so that its local indexing can be rebuilt when it is
     * allocated again.
     *
     * @pre Component and node connections may use subsystem-local indices.
     *
     * @post Component and node connections use their original global system indices,
     *       `global_to_local_map_` is empty, and the subsystem is unallocated.
     *
     * @return 0 on success, otherwise an error code returned while restoring the
     *         global connections.
     */
    int restoreGlobalConnection()
    {
      if (int err_code = mapLocalToGlobal())
      {
        return err_code;
      }

      global_to_local_map_.clear();

      allocated_ = false;

      return 0;
    }

    const std::vector<IdxT>& getExternalDataIndices() const
    {
      return external_data_indices_;
    }

    VectorT& getExternalDataY()
    {
      return y_ext_data_;
    }

    VectorT& getExternalDataYP()
    {
      return yp_ext_data_;
    }

    VectorT& getExternalDataF()
    {
      return f_ext_data_;
    }

    void setForcingFunction(TimeFunction function)
    {
      forcing_function_ = std::move(function);
    }

    const std::unordered_map<IdxT, IdxT>& getInternalMap() const
    {
      return global_to_local_map_;
    }

  private:
    /**
     * @brief Restores subsystem-local connection indices to global system indices.
     *
     * Replaces the local connection indices stored by each component and node with
     * their corresponding global system indices. Each local index is mapped back to
     * its original global connection using the subsystem connection table.
     *
     * This reverses the local indexing established during subsystem allocation
     * without changing the underlying component or node connectivity.
     *
     * @pre Component and node connections use subsystem-local indices.
     *
     * @post Component and node connections use their original global system indices.
     *
     * @return 0 on success.
     */
    int mapLocalToGlobal()
    {
      if (!this->isAllocated())
      {
        return 0;
      }

      IdxT local_connection;

      for (component_type* component : components_)
      {
        for (IdxT i = 0; i < component->size(); i++)
        {
          local_connection = component->getNodeConnection(i);

          if (local_connection != INVALID_INDEX<IdxT>)
          {
            component->setConnectionNodes(i, this->getNodeConnection(local_connection));
          }
        }
      }

      for (node_type* node : nodes_)
      {
        for (IdxT i = 0; i < node->size(); i++)
        {
          local_connection = node->getNodeConnection(i).idx_;

          if (local_connection != INVALID_INDEX<IdxT>)
          {
            node->setConnectionNodes(i, this->getNodeConnection(local_connection));
          }
        }
      }

      return 0;
    }

    /**
     * @brief Builds the global-to-local connection map for the subsystem.
     *
     * Assigns a local connection index to every variable referenced by the subsystem.
     * Variables owned by components and nodes in the subsystem form the internal
     * variables of the subsystem and are assigned the first N local indices [0, N).
     * Variables required by the subsystem but owned outside it are then assigned
     * local indices starting at N.
     *
     * A variable that is external to a component may still be internal to the
     * subsystem if it is owned by another component or node in the same subsystem.
     * Such variables retain the local indices assigned during the internal passes.
     *
     * @pre Component and node connections use global system indices.
     *
     * @post `global_to_local_map_` maps each referenced global connection to its
     *       subsystem-local index, with internal variables followed by external
     *       variables.
     *
     * @return 0 on success.
     */
    int buildConnectionMap()
    {

      global_to_local_map_.clear();

      // Collect each component's internal variables and assign them local indices,
      // since they make up the internal variables of the subsystem.
      IdxT local_connection = 0;
      for (component_type* comp : components_)
      {
        const auto& extern_indices = comp->getExternIndices();

        for (IdxT i = 0; i < comp->size(); i++)
        {
          const IdxT global_connection = comp->getNodeConnection(i);

          if (global_connection != INVALID_INDEX<IdxT> && !extern_indices.contains(i))
          {
            global_to_local_map_[global_connection] = local_connection++;
          }
        }
      }

      // Node variables are also owned by the subsystem and therefore complete
      // the set of subsystem-internal variables.
      for (node_type* node : nodes_)
      {
        for (IdxT i = 0; i < node->size(); i++)
        {
          const IdxT global_connection = node->getNodeConnection(i).idx_;

          if (global_connection != INVALID_INDEX<IdxT>)
          {
            global_to_local_map_[global_connection] = local_connection++;
          }
        }
      }

      // Everything assigned so far is internal to the subsystem.
      const IdxT internal_size = local_connection;

      // Finally assign local indices to variables owned outside the subsystem.
      for (component_type* comp : components_)
      {
        auto extern_indices = comp->getExternIndices();

        for (IdxT j = 0; j < comp->size(); j++)
        {
          if (extern_indices.contains(j))
          {
            const IdxT global_connection = comp->getNodeConnection(j);

            if (global_connection != INVALID_INDEX<IdxT> && !global_to_local_map_.contains(global_connection))
            {
              global_to_local_map_[global_connection] = local_connection++;
            }
          }
        }
      }

      n_intern_ = static_cast<size_t>(internal_size);
      n_extern_ = static_cast<size_t>(local_connection - internal_size);
      size_     = local_connection;

      return 0;
    }

    /**
     *@brief Maps global system connection indices to local subsystem connection indices.
     */
    std::unordered_map<IdxT, IdxT> global_to_local_map_;

    /**
     * @brief Global system index corresponding to each entry in the external subsystem vectors.
     */
    std::vector<IdxT> external_data_indices_;

    /**
     * @brief subsystem external state data.
     */
    VectorT y_ext_data_;

    /**
     * @brief subsystem external state derivative
     */
    VectorT yp_ext_data_;

    /**
     * @brief subsystem external residual data
     */
    VectorT f_ext_data_;

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

  }; // class SubsystemModel

} // namespace GridKit