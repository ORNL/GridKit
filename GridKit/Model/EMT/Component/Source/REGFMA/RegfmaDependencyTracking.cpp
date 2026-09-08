#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>

#include "RegfmaImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    int Regfma<scalar_type, index_type>::assembleJacobian(RealT, RealT)
    {
      return 0;
    }

    template class Regfma<DependencyTracking::Variable, long int>;
    template class Regfma<DependencyTracking::Variable, size_t>;
  } // namespace EMT
} // namespace GridKit
