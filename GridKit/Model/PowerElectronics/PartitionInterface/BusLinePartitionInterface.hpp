/**
 * @file BusLinePartitionInterface.hpp
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
  namespace PowerElectronics
  {
    /**
     * @brief Partition interface between a bus and a line.
     *
     * The BusLinePartitionInterface provides a mechanism for a partition to access
     * the current injection it requires from a line located in a neighboring
     * partition, without evaluating the line within the local partition.
     *
     * The interface tracks four variables associated with the bus-line connection.
     * The first two variables belong to the bus participating in the partition
     * boundary, while the last two belong to the line connected to that bus. All
     * four variables are treated as external variables by the interface.
     *
     * During partitioning, the interface is added to the partition containing the
     * bus. Since the first two interface variables belong to components within the
     * partition, the subsystem maps them to its local variables. The last two
     * variables belong to the line in the neighboring partition and are therefore
     * identified by the subsystem as external variables.
     *
     * During subsystem evaluation, the interface uses these external line variables
     * to provide the corresponding current injection to the bus equations. The
     * current contribution is added to or subtracted from the bus residual according
     * to which end of the line is connected to the bus.
     *
     * @tparam ScalarT Scalar type used for state and residual values.
     * @tparam IdxT Index type used for component and variable indices.
     */
    template <class ScalarT, typename IdxT>
    class BusLinePartitionInterface : public PartitionInterface<ScalarT, IdxT>
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
      BusLinePartitionInterface(node_type* bus, component_type* component, IdxT id);
      virtual ~BusLinePartitionInterface();

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
      /**
       * @brief Owned copy of the line component across the partition boundary.
       */
      std::unique_ptr<component_type> component_;

      /**
       * @brief Non-owning pointer to the bus within the local partition.
       */
      node_type* bus_;

      /**
       * @brief True if the bus is connected to the left end of the line.
       */
      bool is_left_of_Line_;
    };
  } // namespace PowerElectronics
} // namespace GridKit
