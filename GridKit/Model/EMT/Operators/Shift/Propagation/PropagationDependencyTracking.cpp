#include "PropagationImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    template class Propagation<DependencyTracking::Variable, long int>;
    template class Propagation<DependencyTracking::Variable, size_t>;
  } // namespace EMT
} // namespace GridKit
