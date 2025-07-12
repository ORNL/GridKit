/**
 * @file SignalNode model implementation.
 */
#include "SignalNodeImpl.hpp"


namespace GridKit
{
  namespace PhasorDynamics
  {
    template class SignalNode<double, size_t>;
    template class SignalNode<double, long>;

  } // namespace PhasorDynamics
} // namespace GridKit
