/**
 * @file NortonDependencyTracking.cpp
 * @brief Template instantiations for the EMT Norton model.
 */

#include "NortonImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    int Norton<scalar_type, index_type>::assembleJacobian(RealT, RealT)
    {
      return 0;
    }

    template class Norton<DependencyTracking::Variable, long int>;
    template class Norton<DependencyTracking::Variable, size_t>;
  } // namespace EMT
} // namespace GridKit
