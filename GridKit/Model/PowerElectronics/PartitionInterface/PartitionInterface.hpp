

#pragma once

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Model/PowerElectronics/CircuitComponent.hpp>

namespace GridKit
{
  /*!
   * @brief Base class for partition interface components.
   *
   * PartitionInterface provides the infrastructure needed by components
   * that exchange information across partition boundaries during
   * co-simulation. It stores the values of external variables received
   * from neighboring partitions together with the corresponding global
   * indices in the original unpartitioned system.
   *
   * Derived classes are responsible for implementing the specific
   * interface behavior associated with buses, circuit components, or
   * other partition boundary entities.
   */
  template <class ScalarT, typename IdxT>
  class PartitionInterface : public CircuitComponent<ScalarT, IdxT>
  {
  public:
    PartitionInterface() = default;

    ~PartitionInterface() = default;

  protected:
    size_t bus_port_i_;
    size_t bus_port_j_;
  };

} // namespace GridKit
