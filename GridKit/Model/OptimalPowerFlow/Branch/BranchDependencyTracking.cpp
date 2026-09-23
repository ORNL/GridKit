/**
 * @file BranchDependencyTracking.cpp
 * @brief Optimal power flow branch instantiation for `DependencyTracking::Variable`.
 */

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>

#include "BranchImpl.hpp"

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    // Available template instantiations
    template class ComponentModel<Branch<DependencyTracking::Variable, size_t>, DependencyTracking::Variable, size_t>;
    template class Branch<DependencyTracking::Variable, size_t>;
  } // namespace OptimalPowerFlow
} // namespace GridKit
