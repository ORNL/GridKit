
#include "BusImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    // Available template instantiations
    template class Bus<DependencyTracking::Variable, long int>;
    template class Bus<DependencyTracking::Variable, size_t>;

  } // namespace EMT
} // namespace GridKit
