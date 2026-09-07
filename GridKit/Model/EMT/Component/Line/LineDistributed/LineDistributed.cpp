#include "LineDistributedImpl.hpp"
#include "LineDistributedJacobian.hpp"

namespace GridKit
{
  namespace EMT
  {
    // Available template instantiations
    template class LineDistributed<double, long int>;
    template class LineDistributed<double, size_t>;
  } // namespace EMT
} // namespace GridKit
