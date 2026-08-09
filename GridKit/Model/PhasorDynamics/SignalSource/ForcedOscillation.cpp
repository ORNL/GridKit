/**
 * @file ForcedOscillation.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Plain-scalar instantiations for the forced-oscillation signal source.
 */

#include "ForcedOscillationImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    template class ForcedOscillation<double, long int>;
    template class ForcedOscillation<double, size_t>;
  } // namespace PhasorDynamics
} // namespace GridKit
