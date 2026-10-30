/**
 * @file GasBusInfiniteDependencyTracking.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 */

#include "GasBusInfiniteImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Gas
    {
      // Available template instantiations
      template class GasBusInfinite<DependencyTracking::Variable, long int>;
      template class GasBusInfinite<DependencyTracking::Variable, size_t>;
    } // namespace Gas
  } // namespace PhasorDynamics
} // namespace GridKit
