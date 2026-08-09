/**
 * @file SignalNodeJunction.cpp
 * @brief Explicit instantiations of the signal-node junction component.
 */
#include <cstddef>

#include "SignalNodeJunctionImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    template class SignalNodeJunction<double, long int>;
    template class SignalNodeJunction<double, std::size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
