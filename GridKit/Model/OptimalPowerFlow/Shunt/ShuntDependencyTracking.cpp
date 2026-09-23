/**
 * @file ShuntDependencyTracking.cpp
 * @brief Optimal power flow shunt instantiation for `DependencyTracking::Variable`.
 */

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>

#include "ShuntImpl.hpp"

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    // Available template instantiations
    template class ComponentModel<Shunt<DependencyTracking::Variable, size_t>, DependencyTracking::Variable, size_t>;
    template class Shunt<DependencyTracking::Variable, size_t>;
  } // namespace OptimalPowerFlow
} // namespace GridKit
