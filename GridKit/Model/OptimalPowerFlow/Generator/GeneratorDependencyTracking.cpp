/**
 * @file GeneratorDependencyTracking.cpp
 * @brief Optimal power flow generator instantiation for `DependencyTracking::Variable`.
 * Its derivatives come from the system model dependencies.
 */

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>

#include "GeneratorImpl.hpp"

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    template <typename scalar_type, typename index_type>
    int Generator<scalar_type, index_type>::evaluateGradient()
    {
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Generator<scalar_type, index_type>::evaluateJacobian()
    {
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Generator<scalar_type, index_type>::evaluateHessian([[maybe_unused]] RealT sigma)
    {
      return 0;
    }

    // Available template instantiations
    template class Generator<DependencyTracking::Variable, size_t>;
  } // namespace OptimalPowerFlow
} // namespace GridKit
