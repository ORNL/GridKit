#include "AngleImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    template class Angle<DependencyTracking::Variable, long int>;
    template class Angle<DependencyTracking::Variable, size_t>;
  } // namespace EMT
} // namespace GridKit
