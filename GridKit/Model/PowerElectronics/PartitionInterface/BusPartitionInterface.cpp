
#include "BusPartitionInterface.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

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
  BusPartitionInterface<ScalarT, IdxT>::BusPartitionInterface(CircuitComponent<ScalarT, IdxT>& component, IdxT bus_i, IdxT bus_j, IdxT id)
    : component_(component)
  {
    size_     = component.size();
    n_intern_ = 0;
    n_extern_ = static_cast<size_t>(component.size());
    idc_      = id;

    for (IdxT i = 0; i < size_; i++)
    {
      extern_indices_.insert(i);
    }

    bus_i_ = bus_i;
    bus_j_ = bus_j;
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
    y_.resize(static_cast<size_t>(size_));
    yp_.resize(static_cast<size_t>(size_));
    f_.resize(static_cast<size_t>(size_));

    std::fill(f_.begin(), f_.end(), 0);

    component_.allocate();

    for (size_t i = 0; i < static_cast<size_t>(component_.size()); i++)
    {
      if (bus_i_ == component_.getNodeConnection(static_cast<IdxT>(i)))
      {
        bus_port_i_ = i;
      }
      else if (bus_j_ == component_.getNodeConnection(static_cast<IdxT>(i)))
      {
        bus_port_j_ = i;
      }
      IdxT global_idx = component_.getNodeConnection(static_cast<IdxT>(i));
      this->setExternalConnectionNodes(static_cast<IdxT>(i), global_idx);
    }

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
  int BusPartitionInterface<ScalarT, IdxT>::evaluateResidual()
  {
    std::copy(y_.begin(), y_.end(), component_.y().begin());
    std::copy(yp_.begin(), yp_.end(), component_.yp().begin());

    component_.evaluateResidual();

    auto& f = component_.getResidual();

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
    jac_.zeroMatrix();

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
