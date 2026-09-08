#include "FilterImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    int Filter<scalar_type, index_type>::assembleJacobian(RealT, RealT)
    {
      throw std::logic_error("Filter requires an Enzyme-enabled build for Jacobian evaluation");
    }

    template class Filter<double, long int>;
    template class Filter<double, size_t>;
  } // namespace EMT
} // namespace GridKit
