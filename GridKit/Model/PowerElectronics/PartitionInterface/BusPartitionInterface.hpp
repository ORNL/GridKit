

#pragma once

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
    using RealT          = typename CircuitComponent<ScalarT, IdxT>::RealT;

    using CircuitComponent<ScalarT, IdxT>::size_;
    using CircuitComponent<ScalarT, IdxT>::nnz_;
    using CircuitComponent<ScalarT, IdxT>::time_;
    using CircuitComponent<ScalarT, IdxT>::alpha_;
    using CircuitComponent<ScalarT, IdxT>::y_ext_;
    using CircuitComponent<ScalarT, IdxT>::y_int_;
    using CircuitComponent<ScalarT, IdxT>::yp_ext_;
    using CircuitComponent<ScalarT, IdxT>::yp_int_;
    using CircuitComponent<ScalarT, IdxT>::tag_;
    using CircuitComponent<ScalarT, IdxT>::abs_tol_;
    using CircuitComponent<ScalarT, IdxT>::f_ext_;
    using CircuitComponent<ScalarT, IdxT>::f_int_;
    using CircuitComponent<ScalarT, IdxT>::g_;
    using CircuitComponent<ScalarT, IdxT>::yB_;
    using CircuitComponent<ScalarT, IdxT>::ypB_;
    using CircuitComponent<ScalarT, IdxT>::fB_;
    using CircuitComponent<ScalarT, IdxT>::gB_;
    using CircuitComponent<ScalarT, IdxT>::param_;
    using CircuitComponent<ScalarT, IdxT>::idc_;

    using CircuitComponent<ScalarT, IdxT>::extern_indices_;
    using CircuitComponent<ScalarT, IdxT>::n_extern_;
    using CircuitComponent<ScalarT, IdxT>::n_intern_;
    using CircuitComponent<ScalarT, IdxT>::connection_nodes_;

    using PartitionInterface<ScalarT, IdxT>::bus_port_i_;
    using PartitionInterface<ScalarT, IdxT>::bus_port_j_;

  public:
    BusPartitionInterface(node_type& bus, component_type& component, IdxT id);
    virtual ~BusPartitionInterface();

    int allocate();
    int initialize();
    int tagDifferentiable();
    int setAbsoluteTolerance(RealT);
    int evaluateInternalResidual() final;
    int evaluateExternalResidual() final;
    int evaluateJacobian();
    int evaluateIntegrand();

    int initializeAdjoint();
    int evaluateAdjointResidual();
    // int evaluateAdjointJacobian();
    int evaluateAdjointIntegrand();
    int updateComponentPointers(ScalarT* residual);

  private:
    component_type& component_;
    node_type&      bus_;

    std::unique_ptr<ScalarT[]> y_ptr;
    std::unique_ptr<ScalarT[]> yp_ptr;
    std::unique_ptr<ScalarT[]> f_ptr;

    ScalarT dummy_residual_{};

    std::vector<IdxT> jac_map_;
  };
} // namespace GridKit