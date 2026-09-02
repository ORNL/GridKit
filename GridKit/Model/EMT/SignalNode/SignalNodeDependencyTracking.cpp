/**
 * @file SignalNode model implementation.
 */
#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>

#include "SignalNodeImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    template class SignalNode<DependencyTracking::Variable, size_t>;
    template class SignalNode<DependencyTracking::Variable, long>;

  } // namespace EMT
} // namespace GridKit
