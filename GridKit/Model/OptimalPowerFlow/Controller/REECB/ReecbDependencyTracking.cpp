/**
 * @file ReecbDependencyTracking.cpp
 * @brief Optimal power flow REECB controller instantiation for `DependencyTracking::Variable`.
 * Its derivatives come from the system model dependencies.
 */

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>

#include "ReecbImpl.hpp"

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    namespace Controller
    {
      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::evaluateGradient()
      {
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::evaluateJacobian()
      {
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::evaluateHessian([[maybe_unused]] RealT sigma)
      {
        return 0;
      }

      // Available template instantiations
      template class Reecb<DependencyTracking::Variable, size_t>;
    } // namespace Controller
  } // namespace OptimalPowerFlow
} // namespace GridKit
