/**
 * @file Tgov1DependencyTracking.cpp
 * @brief Optimal power flow TGOV1 governor instantiation for `DependencyTracking::Variable`.
 * Its derivatives come from the system model dependencies.
 */

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>

#include "Tgov1Impl.hpp"

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    namespace Governor
    {
      template <typename scalar_type, typename index_type>
      int Tgov1<scalar_type, index_type>::evaluateGradient()
      {
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Tgov1<scalar_type, index_type>::evaluateJacobian()
      {
        return 0;
      }

      // Available template instantiations
      template class Tgov1<DependencyTracking::Variable, size_t>;
    } // namespace Governor
  } // namespace OptimalPowerFlow
} // namespace GridKit
