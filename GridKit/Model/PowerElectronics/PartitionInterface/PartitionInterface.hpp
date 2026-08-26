

#pragma once

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Model/PowerElectronics/CircuitComponent.hpp>

namespace GridKit
{
  /*!
   * @brief Base class for partition interface components.
   *
   * A partition interface is simply a component that is designed to aid in
   * marking partition boundary to facilitate partion evaluation and
   * Jacobian assembly. All partition interface that will be implemented
   * in the future should inherit from this class and all common features
   * should be factored here.
   */
  template <class ScalarT, typename IdxT>
  class PartitionInterface : public CircuitComponent<ScalarT, IdxT>
  {
  public:
    PartitionInterface() = default;

    ~PartitionInterface() = default;
  };

} // namespace GridKit
