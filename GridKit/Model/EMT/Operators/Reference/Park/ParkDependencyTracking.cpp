#include "ParkImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    template class Park<DependencyTracking::Variable, long int>;
    template class Park<DependencyTracking::Variable, size_t>;
  } // namespace EMT
} // namespace GridKit
