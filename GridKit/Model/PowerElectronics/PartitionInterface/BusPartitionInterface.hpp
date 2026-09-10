/**
 * @file BusPartitionInterface.hpp
 * @author Abdourahman Barry (abdourahman@vt.edu)
 *
 */

#pragma once

#include <memory>

#include <GridKit/Model/PowerElectronics/CircuitComponent.hpp>
#include <GridKit/Model/PowerElectronics/NodeBase.hpp>
#include <GridKit/Model/PowerElectronics/PartitionInterface/PartitionInterface.hpp>

namespace GridKit
{
  /**
   * @brief Partition interface that provides a wrapped component's current
   *        injection contributions to a bus.
   *
   * A BusPartitionInterface is used when a bus and one of its connected
   * components belong to different partitions. The interface represents the
   * component within the partition containing the bus and provides the current
   * injection contributions that the wrapped component makes to that bus.
   *
   * The interface has the same number of variables as the wrapped component,
   * with a one-to-one correspondence between interface variables and wrapped
   * component variables. However, all variables are external from the
   * interface's point of view, regardless of whether the corresponding variable
   * is internal or external to the wrapped component:
   *
   * @code
   * Interface variables:      [ext_, ext_, ext_, ext_, ext_, ext_, ext_]
   *                              |     |     |     |     |     |     |
   * Wrapped component:        [int_, int_, ext_, int_, ext_, ext_, ext_]
   * @endcode
   *
   * This allows the partition containing the wrapped component to provide all
   * state information required to reconstruct and evaluate the component on the
   * bus side of the partition boundary. Variables that are internal to the
   * wrapped component are stored in private interface storage during evaluation,
   * while variables that are external to the wrapped component are connected
   * directly to the corresponding interface data.
   *
   * After evaluating the wrapped component, the interface extracts the residual
   * and Jacobian contributions associated with the bus and exposes them to the
   * partition containing the bus.
   *
   * * @note The interface residual contributions are zero for all variables except
   *       those connected to the bus. Only variables connected to the bus expose
   *       the wrapped component's contribution to the partition.
   *
   * @note This interface must be added to the partition containing the bus
   *       to ensure that the system is properly partitioned.
   *
   * @tparam ScalarT Scalar type used by the model.
   * @tparam IdxT Index type used for variables and connections.
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

  public:
    BusPartitionInterface(node_type* bus, component_type* component, IdxT id);
    virtual ~BusPartitionInterface();

    int allocate() final;
    int initialize() final;
    int tagDifferentiable() final;
    int setAbsoluteTolerance(RealT) final;
    int evaluateInternalResidual() final;
    int evaluateExternalResidual() final;
    int evaluateJacobian() final;
    int evaluateIntegrand() final;

    int initializeAdjoint() final;
    int evaluateAdjointResidual() final;
    // int evaluateAdjointJacobian();
    int evaluateAdjointIntegrand() final;

  private:
    int updateComponentPointers();

    /**
     * @brief Circuit component participating in the partition split.
     */
    component_type* component_;

    /**
     * @brief Bus associated with the partition boundary.
     */
    node_type* bus_;

    /**
     * @brief Storage for the wrapped component's internal state variables.
     */
    std::unique_ptr<ScalarT[]> component_y_int_;

    /**
     * @brief Storage for the wrapped component's internal state derivatives.
     */
    std::unique_ptr<ScalarT[]> component_yp_int_;

    /**
     * @brief Storage for the wrapped component's internal residual values.
     */
    std::unique_ptr<ScalarT[]> component_f_int_;

    /**
     * @brief Storage for the wrapped component's external residual values.
     */
    std::unique_ptr<ScalarT[]> component_f_ext_;

    /**
     * @brief Wrapped component ports that receive variables from the bus.
     */
    std::vector<size_t> bus_input_ports_;

    /**
     * @brief Wrapped component ports whose residual contributions are exposed
     *        to the bus.
     */
    std::vector<size_t> bus_output_ports_;

    /**
     * @brief Mapping from interface Jacobian entries to wrapped component
     *        Jacobian entries.
     */
    std::vector<IdxT> jac_map_;

    /**
     * @brief Cached lookup indicating whether each wrapped component variable is external.
     *
     * Avoids repeatedly searching `extern_indices_` when determining whether
     * a wrapped component variable is internal or external.
     */
    std::vector<bool> is_external_;
  };
} // namespace GridKit
