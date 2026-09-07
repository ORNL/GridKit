#include "LineDistributedImpl.hpp"
#include "LineDistributedJacobian.hpp"

namespace GridKit
{
  namespace EMT
  {
    // Available template instantiations
    template class LineDistributed<DependencyTracking::Variable, long int>;
    template class LineDistributed<DependencyTracking::Variable, size_t>;
  } // namespace EMT
} // namespace GridKit
