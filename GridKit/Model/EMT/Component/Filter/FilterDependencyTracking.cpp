#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>

#include "FilterImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    int Filter<scalar_type, index_type>::assembleJacobian(RealT, RealT)
    {
      return 0;
    }

    template class Filter<DependencyTracking::Variable, long int>;
    template class Filter<DependencyTracking::Variable, size_t>;
  } // namespace EMT
} // namespace GridKit
