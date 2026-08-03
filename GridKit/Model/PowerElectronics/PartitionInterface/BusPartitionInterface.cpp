
#include "BusPartitionInterface.hpp"

#include <cassert>
#include <cmath>
#include <cstddef>

#include "GridKit/Model/PowerElectronics/ExternalConnection.hpp"

namespace GridKit
{

  /*!
   * @brief Construct a bus-level partition interface.
   *
   * This interface exposes only the two bus variables associated with a
   * component that has been separated from the rest of the system during
   * partitioning. Internally, it maintains a copy of the original component
   * and acts as a proxy between the co-simulation algorithm and the component
   * model.
   *
   * The bus variables become the internal variables of this interface,
   * while all remaining component variables are treated as external data
   * supplied by neighboring partitions.
   *
   * @param component Circuit component associated with this interface.
   * @param id Unique component identifier.
   */
  template <class ScalarT, typename IdxT>
  BusPartitionInterface<ScalarT, IdxT>::BusPartitionInterface(node_type& bus, component_type& component, IdxT id)
    : component_(component),
      bus_(bus)
  {
    size_     = component_.size();
    n_intern_ = 0;
    n_extern_ = static_cast<size_t>(component_.size());
    idc_      = id;

    for (IdxT i = 0; i < size_; i++)
    {
      extern_indices_.insert(i);
    }

    // The subsystem needs access to the interface connection indices before
    // the interface is fully allocated.
    connection_nodes_ = std::make_unique<IdxT[]>(static_cast<size_t>(size_));

    const IdxT bus_i = bus_.getNodeConnection(0).idx_;
    const IdxT bus_j = bus_.getNodeConnection(1).idx_;

    bool port_i_set = false;
    bool port_j_set = false;

    // Copy the wrapped component topology and identify the two bus ports.
    for (size_t i = 0; i < static_cast<size_t>(size_); ++i)
    {
      const IdxT connection_index = component_.getNodeConnection(i);

      this->setConnectionNodes(i, connection_index);

      if (connection_index == bus_i)
      {
        bus_port_i_ = i;
        port_i_set  = true;
      }
      else if (connection_index == bus_j)
      {
        bus_port_j_ = i;
        port_j_set  = true;
      }
    }

    if (!port_i_set || !port_j_set)
    {
      std::cerr << "ERROR: Invalid partition interface detected. "
                << "Bus(ID=" << bus_.busID()
                << "), Component(ID=" << component_.getIDcomponent()
                << "). Please verify connection-node mappings and "
                   "internal/external index assignments."
                << std::endl;
    }

    const IdxT* cooRow = component_.jacobianCooRows();
    const IdxT* cooCol = component_.jacobianCooCols();

    nnz_ = 0;
    for (IdxT k = 0; k < component_.nnz(); ++k)
    {
      const IdxT row_node = component_.getNodeConnection(static_cast<size_t>(cooRow[k]));
      const IdxT col_node = component_.getNodeConnection(static_cast<size_t>(cooCol[k]));

      const bool row_is_bus = (row_node == bus_i || row_node == bus_j);
      const bool col_is_bus = (col_node == bus_i || col_node == bus_j);

      if (row_is_bus && col_is_bus)
      {
        ++nnz_;
        jac_map_.push_back(k);
      }
    }
  }

  template <class ScalarT, typename IdxT>
  BusPartitionInterface<ScalarT, IdxT>::~BusPartitionInterface()
  {
  }

  /*!
   * @brief allocate method
   */
  template <class ScalarT, typename IdxT>
  int BusPartitionInterface<ScalarT, IdxT>::allocate()
  {

    CircuitComponent<ScalarT, IdxT>::allocate();

    const auto n = component_.getInternalSize();

    y_ptr  = std::make_unique<ScalarT[]>(n);
    yp_ptr = std::make_unique<ScalarT[]>(n);
    f_ptr  = std::make_unique<ScalarT[]>(n);

    component_.setInternalPointer(y_ptr.get());
    component_.setInternalDerivativePointer(yp_ptr.get());
    component_.setInternalResidualPointer(f_ptr.get());

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
    std::vector<ScalarT> component_residual(static_cast<size_t>(component_.size()));

    updateComponentPointers(component_residual.data());

    if (int err_code = component_.evaluateResidual())
    {
      return err_code;
    }

    *f_ext_[bus_port_i_] += component_residual[bus_port_i_];
    *f_ext_[bus_port_j_] += component_residual[bus_port_j_];

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

    component_.evaluateJacobian();

    const IdxT*  cooRows = component_.jacobianCooRows();
    const IdxT*  cooCols = component_.jacobianCooCols();
    const RealT* cooVals = component_.jacobianCooValues();

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

    const auto& external_indices = component_.getExternIndices();

    for (size_t i = 0; i < static_cast<size_t>(component_.size()); ++i)
    {
      if (external_indices.contains(static_cast<IdxT>(i)))
      {
        ExternalConnection<ScalarT, IdxT> connection{
            .y_   = y_ext_[i],
            .yp_  = yp_ext_[i],
            .f_   = residual ? &residual[i] : &dummy_residual_,
            .idx_ = component_.getNodeConnection(i)};

        component_.setExternalConnectionNodes(i, connection);
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