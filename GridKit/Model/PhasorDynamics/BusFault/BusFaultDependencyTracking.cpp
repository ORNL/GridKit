/**
 * @file BusFaultDependencyTracking.cpp
 * @author Slaven Peles (peless@ornl.gov)
 *
 */

#include "BusFaultImpl.hpp"

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
     *
     * @return int - error code, 0 = success
     */
    template <class scalar_type, typename index_type>
    int BusFault<scalar_type, index_type>::evaluateJacobian()
    {
      Log::misc() << "Evaluate DependencyTracking Jacobian for BusFault...\n";
      Log::misc() << "Jacobian evaluation is experimental!\n";

      this->constructCsr();

      return 0;
    }

    // Available template instantiations
    template class BusFault<DependencyTracking::Variable, long int>;
    template class BusFault<DependencyTracking::Variable, size_t>;
  } // namespace PhasorDynamics
} // namespace GridKit
