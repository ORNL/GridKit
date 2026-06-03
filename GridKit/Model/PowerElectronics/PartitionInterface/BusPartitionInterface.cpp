
#include "BusPartitionInterface.hpp"

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
    size_           = 2;
    n_intern_       = 0;
    n_extern_       = 2;
    extern_indices_ = {0, 1};
    idc_            = id;

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

    size_t interface_extern_size = static_cast<size_t>(component_.size() - size_);
    interface_partition_externals_.resize(interface_extern_size);
    external_data_y_.resize(interface_extern_size, 0);
    external_data_yp_.resize(interface_extern_size, 0);

    component_.allocate();

    size_t counter = 0;
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
      else
      {
        interface_partition_externals_[counter] = component_.getNodeConnection(static_cast<IdxT>(i));
        counter++;
      }
    }

    this->setExternalConnectionNodes(0, bus_i_);
    this->setExternalConnectionNodes(1, bus_j_);

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
    auto& y  = component_.y();
    auto& yp = component_.yp();

    y[bus_port_i_] = y_[0];
    y[bus_port_j_] = y_[1];

    yp[bus_port_i_] = yp_[0];
    yp[bus_port_j_] = yp_[1];

    size_t counter = 0;
    for (size_t i = 0; i < static_cast<size_t>(component_.size()); i++)
    {
      if (bus_port_i_ != i && bus_port_j_ != i)
      {
        y[i] = external_data_y_[counter];
        counter++;
      }
    }

    component_.evaluateResidual();

    auto& f = component_.getResidual();

    f_[0] = f[bus_port_i_];
    f_[1] = f[bus_port_j_];

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
