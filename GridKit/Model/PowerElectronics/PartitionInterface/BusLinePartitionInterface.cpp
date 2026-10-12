
#include "BusLinePartitionInterface.hpp"

#include <cassert>
#include <cmath>
#include <cstddef>

#include <GridKit/Model/PowerElectronics/Bus/MicrogridBus.hpp>
#include <GridKit/Model/PowerElectronics/ExternalConnection.hpp>
#include <GridKit/Model/PowerElectronics/MicrogridLine/MicrogridLine.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

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
     * @param bus Bus associated with the partition interface.
     * @param component Line component across the partition boundary.
     * @param id Unique identifier for the partition interface.
     */
    template <class ScalarT, typename IdxT>
    BusLinePartitionInterface<ScalarT, IdxT>::BusLinePartitionInterface(node_type* bus, component_type* component, IdxT id)
      : component_(component->clone()),
        bus_(bus)
    {
      if (dynamic_cast<MicrogridBus<ScalarT, IdxT>*>(bus) == nullptr)
      {
        throw std::invalid_argument("BusLinePartitionInterface requires a MicrogridBus.");
      }

      if (dynamic_cast<MicrogridLine<ScalarT, IdxT>*>(component) == nullptr)
      {
        throw std::invalid_argument("BusLinePartitionInterface requires a MicrogridLine.");
      }

      size_           = 4;
      n_intern_       = 0;
      n_extern_       = 4;
      idc_            = id;
      extern_indices_ = {0, 1, 2, 3};
      nnz_            = 2;
    }

    template <class ScalarT, typename IdxT>
    BusLinePartitionInterface<ScalarT, IdxT>::~BusLinePartitionInterface()
    {
    }

    template <class ScalarT, typename IdxT>
    int BusLinePartitionInterface<ScalarT, IdxT>::allocate()
    {

      CircuitComponent<ScalarT, IdxT>::allocate();

      // Connection indices associated with the two ends of the line.
      auto comp_bus1_idx_i = component_->getNodeConnection(1);
      auto comp_bus1_idx_j = component_->getNodeConnection(2);
      auto comp_bus2_idx_i = component_->getNodeConnection(3);
      auto comp_bus2_idx_j = component_->getNodeConnection(4);

      // Global indices associated with the local bus.
      auto bus_idx_i = bus_->getNodeConnection(0).idx_;
      auto bus_idx_j = bus_->getNodeConnection(1).idx_;

      // Determine whether the local bus is connected to the left or right end of the line.
      bool is_left_connection  = (comp_bus1_idx_i == bus_idx_i) && (comp_bus1_idx_j == bus_idx_j);
      bool is_right_connection = (comp_bus2_idx_i == bus_idx_i) && (comp_bus2_idx_j == bus_idx_j);

      // Neither end or both ends matched -> invalid interface. Exactly one end should match.
      if (is_left_connection == is_right_connection)
      {
        GridKit::Utilities::Logger::error() << "ERROR: Invalid partition interface detected. "
                                            << "Bus(ID=" << bus_->busID()
                                            << "), Component(ID=" << component_->getIDcomponent()
                                            << "). Malformed topology or no partitionable connection was found."
                                            << std::endl;

        return 1;
      }

      if (is_left_connection)
      {
        // Use the connection indices associated with the first end of the line.
        this->setConnectionNodes(0, component_->getNodeConnection(1));
        this->setConnectionNodes(1, component_->getNodeConnection(2));
      }
      else
      {
        // Use the connection indices associated with the second end of the line.
        this->setConnectionNodes(0, component_->getNodeConnection(3));
        this->setConnectionNodes(1, component_->getNodeConnection(4));
      }

      // Map the line current variables that provide the current injection to the bus.
      this->setConnectionNodes(2, component_->getNodeConnection(5));
      this->setConnectionNodes(3, component_->getNodeConnection(6));

      // Store the partition orientation for residual and Jacobian evaluation.
      is_left_of_Line_ = is_left_connection;

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

      const int sign = is_left_of_Line_ ? -1 : 1;

      *f_ext_[0] += sign * *y_ext_[2];
      *f_ext_[1] += sign * *y_ext_[3];

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

      const int sign = is_left_of_Line_ ? -1 : 1;

      // Jacobian entries corresponding to the two transferred line contributions.
      std::vector<IdxT>  rtemp{0, 1};
      std::vector<IdxT>  ctemp{2, 3};
      std::vector<RealT> valtemp{sign * 1.0, sign * 1.0};

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
