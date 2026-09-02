/**
 * @file SignalNode model implementation.
 */
#include "SignalNodeImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    template class SignalNode<double, size_t>;
    template class SignalNode<double, long>;

  } // namespace EMT
} // namespace GridKit
