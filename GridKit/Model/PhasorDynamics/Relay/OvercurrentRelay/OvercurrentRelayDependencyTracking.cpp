/**
 * @file OvercurrentRelayDependencyTracking.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Dependency-tracking instantiations for the overcurrent relay model.
 */

#include "OvercurrentRelayImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Relay
    {
      template <typename scalar_type, typename index_type>
      int OvercurrentRelay<scalar_type, index_type>::evaluateJacobian()
      {
        this->constructCsr();
        return 0;
      }

      template class OvercurrentRelay<DependencyTracking::Variable, long int>;
      template class OvercurrentRelay<DependencyTracking::Variable, size_t>;
    } // namespace Relay
  } // namespace PhasorDynamics
} // namespace GridKit
