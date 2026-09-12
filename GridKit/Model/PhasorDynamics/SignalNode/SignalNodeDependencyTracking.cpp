/**
 * @file SignalNodeDependencyTracking.cpp
 */

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>

#include "SignalNodeImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    // Available template instantiations
    template class SignalNode<DependencyTracking::Variable, size_t>;
    template class SignalNode<DependencyTracking::Variable, long>;

  } // namespace PhasorDynamics
} // namespace GridKit
