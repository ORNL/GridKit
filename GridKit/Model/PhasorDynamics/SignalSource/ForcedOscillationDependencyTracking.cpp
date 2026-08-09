/**
 * @file ForcedOscillationDependencyTracking.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Dependency-tracking instantiations for the forced-oscillation signal source.
 */

#include "ForcedOscillationImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    template class ForcedOscillation<DependencyTracking::Variable, long int>;
    template class ForcedOscillation<DependencyTracking::Variable, size_t>;
  } // namespace PhasorDynamics
} // namespace GridKit
