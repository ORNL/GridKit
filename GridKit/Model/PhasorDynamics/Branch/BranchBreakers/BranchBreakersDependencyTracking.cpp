/**
 * @file BranchBreakersDependencyTracking.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Dependency-tracking instantiations for the breaker-terminated branch model.
 */

#include "BranchBreakersImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename scalar_type, typename index_type>
    int BranchBreakers<scalar_type, index_type>::evaluateJacobian()
    {
      this->constructCsr();
      return 0;
    }

    template class BranchBreakers<DependencyTracking::Variable, long int>;
    template class BranchBreakers<DependencyTracking::Variable, size_t>;
  } // namespace PhasorDynamics
} // namespace GridKit
