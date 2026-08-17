
#include "BusPartitionInterface.hpp"

#include <cassert>
#include <cmath>
#include <cstddef>

#include <GridKit/Model/PowerElectronics/ExternalConnection.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{

  /**
   * @brief Construct a partition interface between a bus and a circuit component.
   *
   * The interface wraps a copy of a component located across a partition
   * boundary and evaluates the contributions that the component makes to the
   * connected bus. All interface variables are treated as external variables
   * and have the same size as the wrapped component.
   *
   *
   * @param bus Bus associated with the partition interface.
   * @param component Copy of the component across the partition boundary.
   * @param id Unique identifier for the interface.
   */
  template <class ScalarT, typename IdxT>
  BusPartitionInterface<ScalarT, IdxT>::BusPartitionInterface(node_type* bus, component_type* component, IdxT id)
    : component_(component->clone()),
      bus_(bus)
  {
    size_     = component_->size();
    n_intern_ = 0;
    n_extern_ = static_cast<size_t>(component_->size());
    idc_      = id;

    // All variables of the bus interface are external to the interface.
    for (IdxT i = 0; i < size_; i++)
    {
      extern_indices_.insert(i);
    }

    // Map each global bus connection index to its position in the bus.
    std::unordered_map<IdxT, size_t> bus_connections;

    for (size_t i = 0; i < static_cast<size_t>(bus_->size()); ++i)
    {
      bus_connections[bus_->getNodeConnection(i).idx_] = i;
    }

    // The interface Jacobian contains only rows corresponding to variables
    // owned by the bus.
    const IdxT* coo_rows = component_->jacobianCooRows();

    nnz_ = 0;

    for (IdxT k = 0; k < component_->nnz(); ++k)
    {
      const IdxT row_node = component_->getNodeConnection(static_cast<size_t>(coo_rows[k]));

      if (bus_connections.contains(row_node))
      {
        ++nnz_;
        jac_map_.push_back(k); // Keep track of the entries so they can be easily extracted
      }
    }
  }

  template <class ScalarT, typename IdxT>
  BusPartitionInterface<ScalarT, IdxT>::~BusPartitionInterface()
  {
    delete component_;
  }

  /**
   * @brief Allocate storage and initialize the interface mappings.
   *
   * Identifies the external variables of the wrapped component that are
   * connected to the bus and builds the mappings used to transfer residual
   * contributions between the component and the interface. Private storage is
   * also allocated for the internal variables of the wrapped component.
   *
   * @return 0 on success and a nonzero value if the component does not contain
   *         all variables required by the bus interface.
   */
  template <class ScalarT, typename IdxT>
  int BusPartitionInterface<ScalarT, IdxT>::allocate()
  {

    CircuitComponent<ScalarT, IdxT>::allocate();

    // Build a lookup of the global indices belonging to the bus.
    std::unordered_map<IdxT, size_t> bus_connections;

    for (size_t i = 0; i < static_cast<size_t>(bus_->size()); ++i)
    {
      bus_connections[bus_->getNodeConnection(i).idx_] = i;
    }

    const auto& external_indices = component_->getExternIndices();

    is_external_.assign(static_cast<size_t>(size_), false);

    bus_input_ports_.clear();
    bus_output_ports_.clear();

    size_t external_index = 0;

    // Identify which component variables are external and which of those external
    // variables are connected to the bus. bus_input_ports_ stores the corresponding
    // variable position in this interface, while bus_output_ports_ stores its position
    // in the wrapped component's external residual vector.
    for (size_t i = 0; i < static_cast<size_t>(size_); ++i)
    {
      const IdxT connection_index = component_->getNodeConnection(i);

      this->setConnectionNodes(i, connection_index);

      if (!external_indices.contains(static_cast<IdxT>(i)))
      {
        continue;
      }

      is_external_[i] = true;

      if (bus_connections.contains(connection_index))
      {
        bus_input_ports_.push_back(i);
        bus_output_ports_.push_back(external_index);
      }

      ++external_index;
    }

    // A valid bus interface must contain every variable belonging to the bus.
    // If fewer bus variables are found, the wrapped component does not provide
    // the complete coupling required by this interface.
    const size_t bus_size = static_cast<size_t>(bus_->size());

    if (bus_input_ports_.size() != bus_size)
    {
      GridKit::Utilities::Logger::error() << "ERROR: Invalid partition interface detected. "
                                          << "Bus(ID=" << bus_->busID()
                                          << "), Component(ID=" << component_->getIDcomponent()
                                          << "). Expected " << bus_size
                                          << " bus connections, but found "
                                          << bus_input_ports_.size() << "."
                                          << std::endl;

      return 1;
    }

    // The wrapped component is evaluated independently by the interface.
    // Allocate storage for its internal variables and residuals and redirect
    // its internal pointers to this private storage.
    const size_t internal_size = static_cast<size_t>(component_->getInternalSize());

    component_y_int_  = std::make_unique<ScalarT[]>(internal_size);
    component_yp_int_ = std::make_unique<ScalarT[]>(internal_size);
    component_f_int_  = std::make_unique<ScalarT[]>(internal_size);

    component_->setInternalPointer(component_y_int_.get());
    component_->setInternalDerivativePointer(component_yp_int_.get());
    component_->setInternalResidualPointer(component_f_int_.get());

    component_f_ext_ = std::make_unique<ScalarT[]>(component_->getExternSize());

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int BusPartitionInterface<ScalarT, IdxT>::setAbsoluteTolerance(RealT rel_tol)
  {
    abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
    return 0;
  }

  /**
   * @brief Initialize the partition interface.
   *
   * @return 0 on success.
   */
  template <class ScalarT, typename IdxT>
  int BusPartitionInterface<ScalarT, IdxT>::initialize()
  {
    return 0;
  }

  /**
   * @brief Identify differential variables
   */
  template <class ScalarT, typename IdxT>
  int BusPartitionInterface<ScalarT, IdxT>::tagDifferentiable()
  {
    return 0;
  }

  /**
   * @brief Eval Internal Residual
   */
  template <class ScalarT, typename IdxT>
  int BusPartitionInterface<ScalarT, IdxT>::evaluateInternalResidual()
  {
    return 0;
  }

  /**
   * @brief Evaluate the wrapped component's contributions to the bus residual.
   *
   * The wrapped component is evaluated using state information supplied through
   * the interface. Residual contributions associated with bus variables are
   * extracted from the component's external residual and accumulated into the
   * corresponding interface residual entries.
   *
   * @return 0 on success, or the error code returned by the wrapped component.
   */
  template <class ScalarT, typename IdxT>
  int BusPartitionInterface<ScalarT, IdxT>::evaluateExternalResidual()
  {
    std::fill_n(component_f_ext_.get(), component_->getExternSize(), ScalarT{});

    updateComponentPointers();

    if (int err_code = component_->evaluateExternalResidual())
    {
      return err_code;
    }

    // Only contributions associated with the bus are accumulated
    // into the interface residual below.
    for (size_t i = 0; i < bus_input_ports_.size(); ++i)
    {
      *f_ext_[bus_input_ports_[i]] += component_f_ext_[bus_output_ports_[i]];
    }

    return 0;
  }

  /**
   * @brief Evaluate the Jacobian contributions associated with the bus.
   *
   * Evaluates the Jacobian of the wrapped component and extracts the entries
   * whose residual rows correspond to variables belonging to the connected bus.
   * The selected entries are then used to assemble the interface Jacobian.
   *
   * @return 0 on success.
   */
  template <class ScalarT, typename IdxT>
  int BusPartitionInterface<ScalarT, IdxT>::evaluateJacobian()
  {

    this->zeroJacMatrix();

    updateComponentPointers();

    component_->evaluateJacobian();

    const IdxT*  cooRows = component_->jacobianCooRows();
    const IdxT*  cooCols = component_->jacobianCooCols();
    const RealT* cooVals = component_->jacobianCooValues();

    std::vector<IdxT>  rows;
    std::vector<IdxT>  cols;
    std::vector<RealT> vals;

    rows.reserve(jac_map_.size());
    cols.reserve(jac_map_.size());
    vals.reserve(jac_map_.size());

    // Extract only the Jacobian entries whose residual rows belong to the bus.
    for (const IdxT index : jac_map_)
    {
      rows.push_back(cooRows[index]);
      cols.push_back(cooCols[index]);
      vals.push_back(cooVals[index]);
    }

    this->setJacValues(rows, cols, vals);

    return 0;
  }

  /**
   * @brief Update the wrapped component with state and output residual pointers.
   *
   * Reconstructs the state required to evaluate the wrapped component using
   * data supplied through the interface. External component variables are
   * connected directly to the interface data, while internal component
   * variables are copied into the interface's private storage, which the wrapped
   * component already keeps track of for its internal state.
   *
   * @pre The interface must be allocated and its external state pointers must
   *      reference valid state and derivative data.
   *
   * @post The wrapped component is configured with the state, derivative, and
   *       residual pointers required for evaluation.
   *
   * @param residual Storage for the component's external residual contributions,
   *                 or nullptr when residual output is not required.
   *
   * @return 0 on success.
   */
  template <class ScalarT, typename IdxT>
  int BusPartitionInterface<ScalarT, IdxT>::updateComponentPointers()
  {
    size_t internal_index = 0;
    size_t external_index = 0;

    // Reconstruct the state expected by the wrapped component from the
    // interface's external data. Route the external residual to component_f_ext_.
    for (size_t i = 0; i < static_cast<size_t>(component_->size()); ++i)
    {
      if (is_external_[i])
      {
        ExternalConnection<ScalarT, IdxT> connection{
            .y_   = y_ext_[i],
            .yp_  = yp_ext_[i],
            .f_   = &component_f_ext_[external_index],
            .idx_ = component_->getNodeConnection(i)};

        component_->setExternalConnectionNodes(i, connection);
        ++external_index;
      }
      else
      {
        component_y_int_[internal_index]  = *y_ext_[i];
        component_yp_int_[internal_index] = *yp_ext_[i];
        ++internal_index;
      }
    }

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int BusPartitionInterface<ScalarT, IdxT>::evaluateIntegrand()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int BusPartitionInterface<ScalarT, IdxT>::initializeAdjoint()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int BusPartitionInterface<ScalarT, IdxT>::evaluateAdjointResidual()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int BusPartitionInterface<ScalarT, IdxT>::evaluateAdjointIntegrand()
  {
    return 0;
  }

  // Available template instantiations
  template class BusPartitionInterface<double, long int>;
  template class BusPartitionInterface<double, size_t>;
  template class BusPartitionInterface<DependencyTracking::Variable, long int>;
  template class BusPartitionInterface<DependencyTracking::Variable, size_t>;

} // namespace GridKit
