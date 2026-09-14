/**
 * @file BranchDependencyTracking.cpp
 * @author Slaven Peles (peless@ornl.gov)
 *
 */

#include "BranchImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Evaluate DependencyTracking::Variable Jacobian.
     *
     * @note Currently only used for testing.
     *
     * No-op, because branch currently does not own residual equations.
     *
     * @return int - error code, 0 = success
     */
    template <typename scalar_type, typename index_type>
    int Branch<scalar_type, index_type>::evaluateJacobian()
    {
      Log::misc() << "Evaluate DependencyTracking Jacobian for Branch...\n";
      Log::misc() << "Jacobian evaluation is experimental!\n";

      return 0;
    }

    // Available template instantiations
    template class Branch<DependencyTracking::Variable, long int>;
    template class Branch<DependencyTracking::Variable, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
