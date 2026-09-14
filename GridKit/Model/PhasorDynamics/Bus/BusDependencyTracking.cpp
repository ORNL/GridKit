/**
 * @file BusDependencyTracking.cpp
 * @author Slaven Peles (peless@ornl.gov)
 *
 */

#include "BusImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Evaluate DependencyTracking::Variable Jacobian.
     *
     * @note Currently only used for testing.
     *
     * DependencyTracking::Variable stores the Jacobian as dependency maps,
     * updated during calls to evaluateResidual().
     * Not yet implemented for bus residuals.
     *
     * @return int - error code
     */
    template <typename scalar_type, typename index_type>
    int Bus<scalar_type, index_type>::evaluateJacobian()
    {
      Log::misc() << "Evaluate DependencyTracking Jacobian for Bus...\n";
      Log::misc() << "Jacobian evaluation is not implemented!\n";

      return 0;
    }

    // Available template instantiations
    template class BusBase<DependencyTracking::Variable, long int>;
    template class BusBase<DependencyTracking::Variable, size_t>;
    template class Bus<DependencyTracking::Variable, long int>;
    template class Bus<DependencyTracking::Variable, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
