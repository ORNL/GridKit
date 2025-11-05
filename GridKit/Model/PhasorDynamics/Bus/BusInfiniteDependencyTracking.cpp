
#include "BusInfiniteImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    // Available template instantiations
    template class BusInfinite<DependencyTracking::Variable, long int>;
    template class BusInfinite<DependencyTracking::Variable, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
