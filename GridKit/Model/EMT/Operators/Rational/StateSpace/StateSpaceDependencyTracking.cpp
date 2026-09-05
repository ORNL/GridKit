#include "StateSpace.hpp"

namespace GridKit
{
  namespace EMT
  {
    template class StateSpace<DependencyTracking::Variable, long int>;
    template class StateSpace<DependencyTracking::Variable, size_t>;
  } // namespace EMT
} // namespace GridKit
