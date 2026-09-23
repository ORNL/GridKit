/**
 * @file GeneratorDependencyTracking.cpp
 * @brief Optimal power flow generator instantiation for `DependencyTracking::Variable`.
 */

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>

#include "GeneratorImpl.hpp"

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    // Available template instantiations
    template class ComponentModel<Generator<DependencyTracking::Variable, size_t>, DependencyTracking::Variable, size_t>;
    template class Generator<DependencyTracking::Variable, size_t>;
  } // namespace OptimalPowerFlow
} // namespace GridKit
