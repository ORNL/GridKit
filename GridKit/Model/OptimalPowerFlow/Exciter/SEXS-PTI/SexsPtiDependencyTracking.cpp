/**
 * @file SexsPtiDependencyTracking.cpp
 * @brief Optimal power flow SEXS-PTI exciter instantiation for `DependencyTracking::Variable`.
 * Its derivatives come from the system model dependencies.
 */

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>

#include "SexsPtiImpl.hpp"

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    namespace Exciter
    {
      template <typename scalar_type, typename index_type>
      int SexsPti<scalar_type, index_type>::evaluateGradient()
      {
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int SexsPti<scalar_type, index_type>::evaluateJacobian()
      {
        return 0;
      }

      // Available template instantiations
      template class SexsPti<DependencyTracking::Variable, size_t>;
    } // namespace Exciter
  } // namespace OptimalPowerFlow
} // namespace GridKit
