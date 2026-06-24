

#pragma once

#include <cstddef>

#include <GridKit/Model/PowerElectronics/CircuitComponent.hpp>
#include <GridKit/Model/PowerElectronics/NodeBase.hpp>

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

    using component_type = CircuitComponent<ScalarT, IdxT>;
    using node_type      = typename PowerElectronics::NodeBase<ScalarT, IdxT>;

    using CircuitComponent<ScalarT, IdxT>::size_;
    using CircuitComponent<ScalarT, IdxT>::nnz_;
    using CircuitComponent<ScalarT, IdxT>::time_;
    using CircuitComponent<ScalarT, IdxT>::alpha_;
    using CircuitComponent<ScalarT, IdxT>::y_;
    using CircuitComponent<ScalarT, IdxT>::y_int_;
    using CircuitComponent<ScalarT, IdxT>::yp_;
    using CircuitComponent<ScalarT, IdxT>::yp_int_;
    using CircuitComponent<ScalarT, IdxT>::tag_;
    using CircuitComponent<ScalarT, IdxT>::f_;
    using CircuitComponent<ScalarT, IdxT>::f_int_;
    using CircuitComponent<ScalarT, IdxT>::g_;
    using CircuitComponent<ScalarT, IdxT>::yB_;
    using CircuitComponent<ScalarT, IdxT>::ypB_;
    using CircuitComponent<ScalarT, IdxT>::fB_;
    using CircuitComponent<ScalarT, IdxT>::gB_;
    using CircuitComponent<ScalarT, IdxT>::param_;
    using CircuitComponent<ScalarT, IdxT>::idc_;

    using PartitionInterface<ScalarT, IdxT>::bus_port_i_;
    using PartitionInterface<ScalarT, IdxT>::bus_port_j_;

    using CircuitComponent<ScalarT, IdxT>::extern_indices_;
    using CircuitComponent<ScalarT, IdxT>::n_extern_;
    using CircuitComponent<ScalarT, IdxT>::n_intern_;

  public:
    BusPartitionInterface(node_type& bus, component_type& component, IdxT id);
    virtual ~BusPartitionInterface();

    int allocate();
    int initialize();
    int tagDifferentiable();
    int evaluateInternalResidual() final;
    int evaluateExternalResidual() final;
    int evaluateJacobian();
    int evaluateIntegrand();

    int initializeAdjoint();
    int evaluateAdjointResidual();
    // int evaluateAdjointJacobian();
    int evaluateAdjointIntegrand();

  private:
    component_type& component_;
    node_type&      bus_;
    ScalarT*        y_ptr;
    ScalarT*        yp_ptr;
    ScalarT*        f_ptr;
  };
} // namespace GridKit
