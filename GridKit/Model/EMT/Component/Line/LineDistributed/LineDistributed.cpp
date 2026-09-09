#include "LineDistributedImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    int LineDistributed<scalar_type, index_type>::assembleJacobian(RealT, RealT)
    {
      throw std::logic_error("LineDistributed Jacobian requires Enzyme");
    }

    // Available template instantiations
    template class LineDistributed<double, long int>;
    template class LineDistributed<double, size_t>;
  } // namespace EMT
} // namespace GridKit
