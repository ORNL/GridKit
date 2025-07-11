#include "BusFaultImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    // Additional template instantiations
    template class BusFault<DependencyTracking::Variable, long int>;
    template class BusFault<DependencyTracking::Variable, size_t>;
  }
}
