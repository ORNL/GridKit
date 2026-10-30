
#include "BusLinePartitionInterface.hpp"

#include <cstddef>
#include <format>

#include <GridKit/Model/PowerElectronics/Bus/MicrogridBus.hpp>
#include <GridKit/Model/PowerElectronics/ExternalConnection.hpp>
#include <GridKit/Model/PowerElectronics/MicrogridLine/MicrogridLine.hpp>

namespace GridKit
{

  namespace PowerElectronics
  {

    /**
     * @brief Construct a partition interface between a bus and a line.
     *
     * The interface consists of four external connections: two bus variables
     * and two line current variables. During partitioning, these connections
     * are mapped to the appropriate subsystem storage.
     *
     * The line currents contribute to the bus residual equations, with signs
     * determined by the orientation of the line relative to the bus.
     *
     * @pre The bus and line connection indices have been allocated.
     *
     * @param bus Bus associated with the partition interface.
     * @param component Line component across the partition boundary.
     * @param id Unique identifier for the partition interface.
     */
    template <class ScalarT, typename IdxT>
    BusLinePartitionInterface<ScalarT, IdxT>::BusLinePartitionInterface(node_type* bus, component_type* component, IdxT id)
    {
      if (dynamic_cast<MicrogridBus<ScalarT, IdxT>*>(bus) == nullptr)
      {
        throw std::invalid_argument("BusLinePartitionInterface: Requires a MicrogridBus.");
      }

      if (dynamic_cast<MicrogridLine<ScalarT, IdxT>*>(component) == nullptr)
      {
        throw std::invalid_argument("BusLinePartitionInterface: Requires a MicrogridLine.");
      }

      if (!bus->isAllocated() or !component->isAllocated())
      {
        throw std::invalid_argument("BusLinePartitionInterface: Malformed network. Ensure components are allocated before partitioning.");
      }

      size_           = 4;
      n_intern_       = 0;
      n_extern_       = 4;
      idc_            = id;
      extern_indices_ = {0, 1, 2, 3};
      nnz_            = 2;

      // Global indices associated with the local bus.
      auto bus_idx_i = bus->getNodeConnection(0).idx_;
      auto bus_idx_j = bus->getNodeConnection(1).idx_;

      const bool is_left_connection  = bus_idx_i == component->getNodeConnection(1) && bus_idx_j == component->getNodeConnection(2);
      const bool is_right_connection = bus_idx_i == component->getNodeConnection(3) && bus_idx_j == component->getNodeConnection(4);

      if (is_left_connection == is_right_connection)
      {
        throw std::invalid_argument(
            std::format(
                "Invalid partition interface detected. Bus(ID={}), Component(ID={}). "
                "Malformed topology or no partitionable connection was found.",
                bus->busID(),
                component->getIDcomponent()));
      }

      sign_ = is_left_connection ? -1.0 : 1.0;
      connection_indices_.resize(static_cast<size_t>(size_));

      connection_indices_[0] = component->getNodeConnection(is_left_connection ? 1 : 3);
      connection_indices_[1] = component->getNodeConnection(is_left_connection ? 2 : 4);
      connection_indices_[2] = component->getNodeConnection(5);
      connection_indices_[3] = component->getNodeConnection(6);
    }

    template <class ScalarT, typename IdxT>
    int BusLinePartitionInterface<ScalarT, IdxT>::allocate()
    {
      Component<ScalarT, IdxT>::allocate();

      this->setConnectionNodes(0, connection_indices_[0]);
      this->setConnectionNodes(1, connection_indices_[1]);
      this->setConnectionNodes(2, connection_indices_[2]);
      this->setConnectionNodes(3, connection_indices_[3]);

      return 0;
    }

    template <class ScalarT, typename IdxT>
    int BusLinePartitionInterface<ScalarT, IdxT>::setAbsoluteTolerance(RealT abs_tol)
    {
      abs_tol_.setToConst(static_cast<ScalarT>(abs_tol));
      return 0;
    }

    /**
     * @brief Initialize the partition interface.
     *
     * @return 0 on success.
     */
    template <class ScalarT, typename IdxT>
    int BusLinePartitionInterface<ScalarT, IdxT>::initialize()
    {
      return 0;
    }

    /**
     * @brief Identify differential variables
     */
    template <class ScalarT, typename IdxT>
    int BusLinePartitionInterface<ScalarT, IdxT>::tagDifferentiable()
    {
      return 0;
    }

    /**
     * @brief Evaluates the internal residual.
     *
     * All variables are external from the interface's point of view, so there
     * are no internal residual contributions to evaluate.
     *
     * @return 0 on success.
     */
    template <class ScalarT, typename IdxT>
    int BusLinePartitionInterface<ScalarT, IdxT>::evaluateInternalResidual()
    {
      return 0;
    }

    /**
     * @brief Evaluate the line contribution to the external bus residual.
     *
     * Adds the two line current residual contributions to the corresponding
     * bus residual equations. The sign of each contribution is determined by
     * whether the bus is connected to the left or right end of the line,
     * consistent with the line current orientation.
     *
     * @return 0 on success.
     */
    template <class ScalarT, typename IdxT>
    int BusLinePartitionInterface<ScalarT, IdxT>::evaluateExternalResidual()
    {
      *f_ext_[0] += sign_ * *y_ext_[2];
      *f_ext_[1] += sign_ * *y_ext_[3];

      return 0;
    }

    /**
     * @brief Evaluate the Jacobian contributions associated with the bus.
     *
     * Evaluates the Jacobian of the entries whose residual rows correspond
     * to variables belonging to the connected bus.
     *
     * @note Some of these entries may not be required to form the Jacobian of the
     *       subsystem, but may be useful for parallel Jacobian assembly.
     *
     * @return 0 on success.
     */
    template <class ScalarT, typename IdxT>
    int BusLinePartitionInterface<ScalarT, IdxT>::evaluateJacobian()
    {
      this->zeroJacMatrix();

      // Jacobian entries corresponding to the two transferred line contributions.
      std::vector<IdxT>  rtemp{0, 1};
      std::vector<IdxT>  ctemp{2, 3};
      std::vector<RealT> valtemp{sign_, sign_};

      this->setJacValues(rtemp, ctemp, valtemp);

      return 0;
    }

    template <class ScalarT, typename IdxT>
    int BusLinePartitionInterface<ScalarT, IdxT>::evaluateIntegrand()
    {
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int BusLinePartitionInterface<ScalarT, IdxT>::initializeAdjoint()
    {
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int BusLinePartitionInterface<ScalarT, IdxT>::evaluateAdjointResidual()
    {
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int BusLinePartitionInterface<ScalarT, IdxT>::evaluateAdjointIntegrand()
    {
      return 0;
    }

    // Available template instantiations
    template class BusLinePartitionInterface<double, long int>;
    template class BusLinePartitionInterface<double, size_t>;
    template class BusLinePartitionInterface<DependencyTracking::Variable, long int>;
    template class BusLinePartitionInterface<DependencyTracking::Variable, size_t>;
  } // namespace PowerElectronics
} // namespace GridKit
