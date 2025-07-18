/**
 * @file SignalNode model implementation.
 */
#include "SignalNodeImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    template class SignalNode<DependencyTracking::Variable, size_t>;
    template class SignalNode<DependencyTracking::Variable, long>;

  } // namespace PhasorDynamics
} // namespace GridKit
