
#include "BusPartitionInterface.hpp"

#include <cassert>
#include <cmath>
#include <cstddef>

#include <GridKit/Model/PowerElectronics/ExternalConnection.hpp>

namespace GridKit
{

  /**
   * @brief Construct a bus partition interface for a circuit component.
   *
   * During partitioning, the interface receives a bus and a component that
   * participate in the partition boundary. The interface is added directly
   * to the partition containing the bus and uses the component model to
   * compute the residual and Jacobian contributions to that bus.
   *
   * The component passed to this constructor must be a copy of the original
   * component, not the original component itself. This allows the interface
   * to modify and evaluate the copied component without affecting the component
   * stored in the neighboring partition.
   *
   * For example:
   * @code
   * CompType* comp_copy = new CompType(*comp);
   * auto* interface = new BusPartitionInterface(bus, comp_copy, id);
   * @endcode
   *
   * @param bus Bus belonging to the partition where the interface is added.
   * @param component Copy of the circuit component across the partition boundary.
   * @param id Unique interface identifier.
   */
  template <class ScalarT, typename IdxT>
  BusPartitionInterface<ScalarT, IdxT>::BusPartitionInterface(node_type* bus, component_type* component, IdxT id)
    : component_(component),
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

    // The subsystem needs access to the interface connection indices before
    // the interface is fully allocated.
    connection_nodes_ = std::make_unique<IdxT[]>(static_cast<size_t>(size_));

    // Global indices of the two variables associated with the bus.
    const IdxT bus_i = bus_->getNodeConnection(0).idx_;
    const IdxT bus_j = bus_->getNodeConnection(1).idx_;

    is_external_.resize(static_cast<size_t>(size_), false);

    bool port_i_set = false;
    bool port_j_set = false;

    const auto& external_indices = component_->getExternIndices();
    IdxT        external_index   = 0;

    // Copy the wrapped component topology and identify the two bus ports and
    // their corresponding output ports in the wrapped component.
    for (size_t i = 0; i < static_cast<size_t>(size_); ++i)
    {
      const IdxT connection_index = component_->getNodeConnection(i);

      this->setConnectionNodes(i, connection_index);

      // Skip indices that are not external to the wrapped component.
      if (!external_indices.contains(static_cast<IdxT>(i)))
      {
        continue;
      }

      // Identify the bus ports and their positions in the component's external output.
      if (connection_index == bus_i)
      {
        bus_port_i_     = i;
        bus_port_out_i_ = static_cast<size_t>(external_index);
        port_i_set      = true;
      }
      else if (connection_index == bus_j)
      {
        bus_port_j_     = i;
        bus_port_out_j_ = static_cast<size_t>(external_index);
        port_j_set      = true;
      }

      ++external_index;
      is_external_[static_cast<size_t>(i)] = true;
    }

    // Report an error if either bus connection is missing from the component.
    if (!port_i_set || !port_j_set)
    {
      std::cerr << "ERROR: Invalid partition interface detected. "
                << "Bus(ID=" << bus_->busID()
                << "), Component(ID=" << component_->getIDcomponent()
                << "). Please verify connection-node mappings and "
                   "internal/external index assignments."
                << std::endl;
    }

    // Record the Jacobian entries whose residual row belongs to either bus
    // variable. jac_map_ stores their positions in the wrapped component's
    // COO Jacobian so that the same entries can be extracted during each Jacobian evaluation.

    const IdxT* cooRow = component_->jacobianCooRows();

    nnz_ = 0;
    for (IdxT k = 0; k < component_->nnz(); ++k)
    {
      const IdxT row_node = component_->getNodeConnection(static_cast<size_t>(cooRow[k]));

      const bool row_is_bus = (row_node == bus_i || row_node == bus_j);

      if (row_is_bus)
      {
        ++nnz_;
        jac_map_.push_back(k);
      }
    }
  }

  template <class ScalarT, typename IdxT>
  BusPartitionInterface<ScalarT, IdxT>::~BusPartitionInterface()
  {
    delete component_;
  }

  /*!
   * @brief allocate method
   */
  template <class ScalarT, typename IdxT>
  int BusPartitionInterface<ScalarT, IdxT>::allocate()
  {

    CircuitComponent<ScalarT, IdxT>::allocate();

    const auto n = component_->getInternalSize();

    y_ptr  = std::make_unique<ScalarT[]>(n);
    yp_ptr = std::make_unique<ScalarT[]>(n);
    f_ptr  = std::make_unique<ScalarT[]>(n);

    component_->setInternalPointer(y_ptr.get());
    component_->setInternalDerivativePointer(yp_ptr.get());
    component_->setInternalResidualPointer(f_ptr.get());

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int BusPartitionInterface<ScalarT, IdxT>::setAbsoluteTolerance(RealT rel_tol)
  {
    abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
    return 0;
  }

  /**
   * Initialization of the grid model
   */
  template <class ScalarT, typename IdxT>
  int BusPartitionInterface<ScalarT, IdxT>::initialize()
  {
    return 0;
  }

  /*
   * \brief Identify differential variables
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

  template <class ScalarT, typename IdxT>
  int BusPartitionInterface<ScalarT, IdxT>::evaluateExternalResidual()
  {
    std::vector<ScalarT> component_ext_residual(static_cast<size_t>(component_->getExternSize()));

    updateComponentPointers(component_ext_residual.data());

    if (int err_code = component_->evaluateResidual())
    {
      return err_code;
    }

    *f_ext_[bus_port_i_] += component_ext_residual[bus_port_out_i_];
    *f_ext_[bus_port_j_] += component_ext_residual[bus_port_out_j_];

    return 0;
  }

  /**
   * @brief Generate Jacobian for
   *
   * @tparam ScalarT
   * @tparam IdxT
   * @return int
   */
  template <class ScalarT, typename IdxT>
  int BusPartitionInterface<ScalarT, IdxT>::evaluateJacobian()
  {

    this->zeroJacMatrix();

    // The Jacobian only requires the component state pointers. Residual
    // contributions are redirected to dummy storage.
    updateComponentPointers(nullptr);

    component_->evaluateJacobian();

    const IdxT*  cooRows = component_->jacobianCooRows();
    const IdxT*  cooCols = component_->jacobianCooCols();
    const RealT* cooVals = component_->jacobianCooValues();

    std::vector<IdxT>  r = {};
    std::vector<IdxT>  c = {};
    std::vector<RealT> v = {};

    for (const auto& index : jac_map_)
    {
      r.push_back(cooRows[index]);
      c.push_back(cooCols[index]);
      v.push_back(cooVals[index]);
    }

    this->setJacValues(r, c, v);

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int BusPartitionInterface<ScalarT, IdxT>::updateComponentPointers(ScalarT* residual)
  {
    size_t internal_index = 0;
    size_t external_index = 0;

    const auto& external_indices = component_->getExternIndices();

    for (size_t i = 0; i < static_cast<size_t>(component_->size()); ++i)
    {
      if (is_external_[i])
      {
        ExternalConnection<ScalarT, IdxT> connection{
            .y_   = y_ext_[i],
            .yp_  = yp_ext_[i],
            .f_   = residual ? &residual[external_index] : &dummy_residual_,
            .idx_ = component_->getNodeConnection(i)};

        component_->setExternalConnectionNodes(i, connection);

        ++external_index;
      }
      else
      {
        y_ptr[internal_index]  = *y_ext_[i];
        yp_ptr[internal_index] = *yp_ext_[i];

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
