/**
 * @file RegcaDependencyTracking.cpp
 * @brief Optimal power flow REGCA converter instantiation for `DependencyTracking::Variable`.
 * Its derivatives come from the system model dependencies.
 */

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>

#include "RegcaImpl.hpp"

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    namespace Converter
    {
      template <typename scalar_type, typename index_type>
      int Regca<scalar_type, index_type>::evaluateGradient()
      {
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Regca<scalar_type, index_type>::evaluateJacobian()
      {
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Regca<scalar_type, index_type>::evaluateHessian([[maybe_unused]] RealT sigma)
      {
        return 0;
      }

      // Available template instantiations
      template class Regca<DependencyTracking::Variable, size_t>;
    } // namespace Converter
  } // namespace OptimalPowerFlow
} // namespace GridKit
