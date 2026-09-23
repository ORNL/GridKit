/**
 * @file BusDependencyTracking.cpp
 * @brief Optimal power flow bus instantiation for `DependencyTracking::Variable`.
 */

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>

#include "BusImpl.hpp"

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    // Available template instantiations
    template class ComponentModel<Bus<DependencyTracking::Variable, size_t>, DependencyTracking::Variable, size_t>;
    template class Bus<DependencyTracking::Variable, size_t>;
  } // namespace OptimalPowerFlow
} // namespace GridKit
