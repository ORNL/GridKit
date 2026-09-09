#include "LineDistributedImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    int LineDistributed<scalar_type, index_type>::assembleJacobian(RealT, RealT)
    {
      return 0;
    }

    // Available template instantiations
    template class LineDistributed<DependencyTracking::Variable, long int>;
    template class LineDistributed<DependencyTracking::Variable, size_t>;
  } // namespace EMT
} // namespace GridKit
