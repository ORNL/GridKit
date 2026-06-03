

#include "ComponentPartitionInterface.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace GridKit
{

  /*!
   * @brief Constructor for a constant MicrogridLoad model
   *
   * Calls default ModelEvaluatorImpl constructor.
   *
   *
   *
   */

  template <class ScalarT, typename IdxT>
  ComponentPartitionInterface<ScalarT, IdxT>::ComponentPartitionInterface(CircuitComponent<ScalarT, IdxT>* component, IdxT bus_i, IdxT bus_j, IdxT id)
    : component_(component)
  {
    // internals [id, iq]
    // externals [\omegaref, vbd_out, vbq_out]
    size_           = component_->size();
    n_intern_       = component_->getInternalSize();
    n_extern_       = component_->getExternSize();
    extern_indices_ = component_->getExternIndices();
    idc_            = id;

    bus_i_ = bus_i;
    bus_j_ = bus_j;
  }

  template <class ScalarT, typename IdxT>
  ComponentPartitionInterface<ScalarT, IdxT>::~ComponentPartitionInterface()
  {
  }

  /*!
   * @brief allocate method computes sparsity pattern of the Jacobian.
   */
  template <class ScalarT, typename IdxT>
  int ComponentPartitionInterface<ScalarT, IdxT>::allocate()
  {

    size_t bus_size = 2;

    y_.resize(static_cast<size_t>(size_));
    yp_.resize(static_cast<size_t>(size_));
    f_.resize(static_cast<size_t>(size_));
    interface_partition_externals_.resize(bus_size);
    external_data_y_.resize(bus_size, 0);
    external_data_yp_.resize(bus_size, 0);

    component_->allocate();

    for (size_t i = 0; i < static_cast<size_t>(size_); i++)
    {
      auto index = component_->getNodeConnection(static_cast<IdxT>(i));
      if (bus_i_ == index)
      {
        bus_port_i_ = i;
      }
      else if (bus_j_ == index)
      {
        bus_port_j_ = i;
      }
      else
      {
        this->setExternalConnectionNodes(static_cast<IdxT>(i), index);
      }
    }

    this->setExternalConnectionNodes(static_cast<IdxT>(bus_port_i_), static_cast<IdxT>(-1));
    this->setExternalConnectionNodes(static_cast<IdxT>(bus_port_j_), static_cast<IdxT>(-1));

    interface_partition_externals_[0] = bus_i_;
    interface_partition_externals_[1] = bus_j_;

    return 0;
  }

  /**
   * Initialization of the grid model
   */
  template <class ScalarT, typename IdxT>
  int ComponentPartitionInterface<ScalarT, IdxT>::initialize()
  {
    return 0;
  }

  /*
   * \brief Identify differential variables
   */
  template <class ScalarT, typename IdxT>
  int ComponentPartitionInterface<ScalarT, IdxT>::tagDifferentiable()
  {
    return 0;
  }

  /**
   * @brief Eval Micro Load
   */
  template <class ScalarT, typename IdxT>
  int ComponentPartitionInterface<ScalarT, IdxT>::evaluateResidual()
  {

    auto& y  = component_->y();
    auto& yp = component_->yp();

    std::copy(y_.begin(), y_.end(), y.begin());
    std::copy(yp_.begin(), yp_.end(), yp.begin());

    y[bus_port_i_] = external_data_y_[0];
    y[bus_port_j_] = external_data_y_[1];

    yp[bus_port_i_] = external_data_yp_[0];
    yp[bus_port_j_] = external_data_yp_[1];

    component_->evaluateResidual();
    auto& f = component_->getResidual();

    std::copy(f.begin(), f.end(), f_.begin());

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
  int ComponentPartitionInterface<ScalarT, IdxT>::evaluateJacobian()
  {
    jac_.zeroMatrix();

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int ComponentPartitionInterface<ScalarT, IdxT>::evaluateIntegrand()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int ComponentPartitionInterface<ScalarT, IdxT>::initializeAdjoint()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int ComponentPartitionInterface<ScalarT, IdxT>::evaluateAdjointResidual()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int ComponentPartitionInterface<ScalarT, IdxT>::evaluateAdjointIntegrand()
  {
    return 0;
  }

  // Available template instantiations
  template class ComponentPartitionInterface<double, long int>;
  template class ComponentPartitionInterface<double, size_t>;
  template class ComponentPartitionInterface<DependencyTracking::Variable, long int>;
  template class ComponentPartitionInterface<DependencyTracking::Variable, size_t>;

} // namespace GridKit
