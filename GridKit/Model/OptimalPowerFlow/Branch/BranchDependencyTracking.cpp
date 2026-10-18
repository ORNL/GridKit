/**
 * @file BranchDependencyTracking.cpp
 * @brief Optimal power flow branch instantiation for `DependencyTracking::Variable`.
 * Its derivatives come from the system model dependencies.
 */

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>

#include "BranchImpl.hpp"

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    template <typename scalar_type, typename index_type>
    int Branch<scalar_type, index_type>::evaluateGradient()
    {
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Branch<scalar_type, index_type>::evaluateJacobian()
    {
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Branch<scalar_type, index_type>::evaluateHessian([[maybe_unused]] RealT sigma)
    {
      return 0;
    }

    // Available template instantiations
    template class Branch<DependencyTracking::Variable, size_t>;
  } // namespace OptimalPowerFlow
} // namespace GridKit
