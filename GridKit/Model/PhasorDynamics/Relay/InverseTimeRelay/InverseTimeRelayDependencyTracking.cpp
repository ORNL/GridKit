/**
 * @file InverseTimeRelayDependencyTracking.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Dependency-tracking instantiations for the InverseTimeRelay model.
 */

#include "InverseTimeRelayImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Relay
    {
      /**
       * @brief Evaluate DependencyTracking::Variable Jacobian.
       *
       * @note Currently only used for testing.
       *
       * DependencyTracking::Variable stores the Jacobian as dependency maps,
       * updated during calls to evaluateResidual().
       */
      template <typename scalar_type, typename index_type>
      int InverseTimeRelay<scalar_type, index_type>::evaluateJacobian()
      {
        this->constructCsr();

        return 0;
      }

      // Available template instantiations
      template class InverseTimeRelay<DependencyTracking::Variable, long int>;
      template class InverseTimeRelay<DependencyTracking::Variable, size_t>;
    } // namespace Relay
  } // namespace PhasorDynamics
} // namespace GridKit
