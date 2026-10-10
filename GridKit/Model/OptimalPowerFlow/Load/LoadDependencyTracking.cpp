/**
 * @file LoadDependencyTracking.cpp
 * @brief Optimal power flow load instantiation for `DependencyTracking::Variable`.
 */

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>

#include "LoadImpl.hpp"

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    // Available template instantiations
    template class ComponentModel<Load<DependencyTracking::Variable, size_t>, DependencyTracking::Variable, size_t>;
    template class Load<DependencyTracking::Variable, size_t>;
  } // namespace OptimalPowerFlow
} // namespace GridKit
