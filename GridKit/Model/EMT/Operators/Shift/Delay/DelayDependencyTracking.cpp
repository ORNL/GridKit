#include "DelayImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    template class Delay<DependencyTracking::Variable, long int>;
    template class Delay<DependencyTracking::Variable, size_t>;
  } // namespace EMT
} // namespace GridKit
