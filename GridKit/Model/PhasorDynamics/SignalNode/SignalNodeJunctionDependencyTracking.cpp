/**
 * @file SignalNodeJunctionDependencyTracking.cpp
 * @brief Dependency-tracking instantiations of the signal-node junction.
 */
#include <cstddef>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>

#include "SignalNodeJunctionImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    template class SignalNodeJunction<DependencyTracking::Variable, long int>;
    template class SignalNodeJunction<DependencyTracking::Variable, std::size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
