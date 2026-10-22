/**
 * @file BranchBreakersDependencyTracking.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Dependency-tracking instantiations for the BranchBreakers model.
 */

#include "BranchBreakersImpl.hpp"

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
     */
    template <typename scalar_type, typename index_type>
    int BranchBreakers<scalar_type, index_type>::evaluateJacobian()
    {
      this->constructCsr();

      return 0;
    }

    // Available template instantiations
    template class BranchBreakers<DependencyTracking::Variable, long int>;
    template class BranchBreakers<DependencyTracking::Variable, size_t>;
  } // namespace PhasorDynamics
} // namespace GridKit
