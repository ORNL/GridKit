
#include "BusPartitionInterface.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>

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

    y_.resize(static_cast<size_t>(size_));
    yp_.resize(static_cast<size_t>(size_));
    f_.resize(static_cast<size_t>(size_));

    std::fill(f_.begin(), f_.end(), 0);
    bool port_i_set = false;
    bool port_j_set = false;

    for (size_t i = 0; i < static_cast<size_t>(size_); i++)
    {
      if (bus_.getNodeConnection(0) == component_.getNodeConnection(static_cast<IdxT>(i)))
      {
        bus_port_i_ = i;
        port_i_set  = true;
      }
      else if (bus_.getNodeConnection(1) == component_.getNodeConnection(static_cast<IdxT>(i)))
      {
        bus_port_j_ = i;
        port_j_set  = true;
      }
      IdxT global_idx = component_.getNodeConnection(static_cast<IdxT>(i));
      this->setExternalConnectionNodes(static_cast<IdxT>(i), global_idx);
    }

    if (!port_i_set || !port_j_set)
    {
      std::cerr << "ERROR: Invalid partition interface detected. "
                << "Bus(ID=" << bus_.busID()
                << "), Component(ID=" << component_.getIDcomponent()
                << "). Please verify connection-node mappings and internal/external index assignments."
                << std::endl;
      assert(false);
    }

    y_ptr  = new ScalarT[component_.getInternalSize()];
    yp_ptr = new ScalarT[component_.getInternalSize()];
    f_ptr  = new ScalarT[component_.getInternalSize()];

    component_.setInternalPointer(y_ptr);
    component_.setInternalDerivativePointer(yp_ptr);
    component_.setInternalResidualPointer(f_ptr);

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
   * @brief Eval Micro Load
   */
  template <class ScalarT, typename IdxT>
  int BusPartitionInterface<ScalarT, IdxT>::evaluateInternalResidual()
  {
    return 0;
  }

  /**
   * @brief Eval Micro Load
   */
  template <class ScalarT, typename IdxT>
  int BusPartitionInterface<ScalarT, IdxT>::evaluateExternalResidual()
  {

    size_t counter_int    = 0;
    size_t counter_ext    = 0;
    auto   extern_indices = component_.getExternIndices();

    for (size_t i = 0; i < static_cast<size_t>(component_.size()); i++)
    {
      if (extern_indices.contains(static_cast<IdxT>(i)))
      {
        component_.y()[counter_ext]  = y_[i];
        component_.yp()[counter_ext] = yp_[i];
        counter_ext++;
      }
      else
      {
        y_ptr[counter_int]  = y_[i];
        yp_ptr[counter_int] = yp_[i];
        counter_int++;
      }
    }

    component_.evaluateExternalResidual();

    auto f = component_.getResidual();

    // TODO: This assumes that external variables are ordered after all internal
    // variables in the local indexing. To make this more robust, we need to get rid of this assumption.
    f_[bus_port_i_] = f[bus_port_i_];
    f_[bus_port_j_] = f[bus_port_j_];
    return 0;
  }

  /**
   * @brief Generate Jacobian for Micro Load
   *
   * @tparam ScalarT
   * @tparam IdxT
   * @return int
   */
  template <class ScalarT, typename IdxT>
  int BusPartitionInterface<ScalarT, IdxT>::evaluateJacobian()
  {
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
