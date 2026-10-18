/**
 * @file LoadDependencyTracking.cpp
 * @brief Optimal power flow load instantiation for `DependencyTracking::Variable`.
 * Its derivatives come from the system model dependencies.
 */

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>

#include "LoadImpl.hpp"

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    template <typename scalar_type, typename index_type>
    int Load<scalar_type, index_type>::evaluateGradient()
    {
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Load<scalar_type, index_type>::evaluateEqualityJacobian()
    {
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Load<scalar_type, index_type>::evaluateInequalityJacobian()
    {
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Load<scalar_type, index_type>::evaluateHessian([[maybe_unused]] RealT        sigma,
                                                       [[maybe_unused]] const RealT* lambda_c,
                                                       [[maybe_unused]] const RealT* lambda_d)
    {
      return 0;
    }

    // Available template instantiations
    template class Load<DependencyTracking::Variable, size_t>;
  } // namespace OptimalPowerFlow
} // namespace GridKit
