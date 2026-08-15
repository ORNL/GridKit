

#pragma once

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Model/PowerElectronics/CircuitComponent.hpp>

namespace GridKit
{
  /*!
   * @brief Base class for partition interface components.
   *
   */
  template <class ScalarT, typename IdxT>
  class PartitionInterface : public CircuitComponent<ScalarT, IdxT>
  {
  public:
    PartitionInterface() = default;

    ~PartitionInterface() = default;
  };

} // namespace GridKit
