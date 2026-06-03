

#pragma once

#include <cstddef>

#include <GridKit/Model/PowerElectronics/CircuitComponent.hpp>

#include "GridKit/Model/PowerElectronics/PartitionInterface/PartitionInterface.hpp"

namespace GridKit
{
  /*!
   * @brief Declaration of a passive BusPartitionInterface class.
   *
   */
  template <class ScalarT, typename IdxT>
  class BusPartitionInterface : public PartitionInterface<ScalarT, IdxT>
  {
    using RealT   = typename CircuitComponent<ScalarT, IdxT>::RealT;
    using MatrixT = typename CircuitComponent<RealT, IdxT>::MatrixT;

    using PartitionInterface<ScalarT, IdxT>::external_data_y_;
    using PartitionInterface<ScalarT, IdxT>::external_data_yp_;
    using PartitionInterface<ScalarT, IdxT>::interface_partition_externals_;
    using PartitionInterface<ScalarT, IdxT>::bus_port_i_;
    using PartitionInterface<ScalarT, IdxT>::bus_port_j_;
    using PartitionInterface<ScalarT, IdxT>::bus_i_;
    using PartitionInterface<ScalarT, IdxT>::bus_j_;

    using CircuitComponent<ScalarT, IdxT>::size_;
    using CircuitComponent<ScalarT, IdxT>::nnz_;
    using CircuitComponent<ScalarT, IdxT>::time_;
    using CircuitComponent<ScalarT, IdxT>::alpha_;
    using CircuitComponent<ScalarT, IdxT>::y_;
    using CircuitComponent<ScalarT, IdxT>::yp_;
    using CircuitComponent<ScalarT, IdxT>::tag_;
    using CircuitComponent<ScalarT, IdxT>::f_;
    using CircuitComponent<ScalarT, IdxT>::g_;
    using CircuitComponent<ScalarT, IdxT>::yB_;
    using CircuitComponent<ScalarT, IdxT>::ypB_;
    using CircuitComponent<ScalarT, IdxT>::fB_;
    using CircuitComponent<ScalarT, IdxT>::gB_;
    using CircuitComponent<ScalarT, IdxT>::jac_;
    using CircuitComponent<ScalarT, IdxT>::param_;
    using CircuitComponent<ScalarT, IdxT>::idc_;

    using CircuitComponent<ScalarT, IdxT>::extern_indices_;
    using CircuitComponent<ScalarT, IdxT>::n_extern_;
    using CircuitComponent<ScalarT, IdxT>::n_intern_;

  public:
    BusPartitionInterface(CircuitComponent<ScalarT, IdxT>& component, IdxT bus_i, IdxT bus_j, IdxT id);
    virtual ~BusPartitionInterface();

    int allocate();
    int initialize();
    int tagDifferentiable();
    int evaluateResidual();
    int evaluateJacobian();
    int evaluateIntegrand();

    int initializeAdjoint();
    int evaluateAdjointResidual();
    // int evaluateAdjointJacobian();
    int evaluateAdjointIntegrand();

  private:
    CircuitComponent<ScalarT, IdxT>& component_;
  };
} // namespace GridKit
