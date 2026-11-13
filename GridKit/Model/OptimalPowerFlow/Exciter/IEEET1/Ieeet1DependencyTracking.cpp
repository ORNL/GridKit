/**
 * @file Ieeet1DependencyTracking.cpp
 * @brief Optimal power flow IEEET1 exciter instantiation for `DependencyTracking::Variable`.
 * Its derivatives come from the system model dependencies.
 */

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>

#include "Ieeet1Impl.hpp"

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    namespace Exciter
    {
      template <typename scalar_type, typename index_type>
      int Ieeet1<scalar_type, index_type>::evaluateGradient()
      {
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Ieeet1<scalar_type, index_type>::evaluateJacobian()
      {
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Ieeet1<scalar_type, index_type>::evaluateHessian([[maybe_unused]] RealT sigma)
      {
        return 0;
      }

      // Available template instantiations
      template class Ieeet1<DependencyTracking::Variable, size_t>;
    } // namespace Exciter
  } // namespace OptimalPowerFlow
} // namespace GridKit
