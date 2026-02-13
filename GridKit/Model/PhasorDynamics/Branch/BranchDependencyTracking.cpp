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
     * @brief Jacobian evaluation not implemented
     *
     * @tparam ScalarT - scalar data type
     * @tparam IdxT    - matrix index data type
     * @return int - error code, 0 = success
     */
    template <class ScalarT, typename IdxT>
    int Branch<ScalarT, IdxT>::evaluateJacobian()
    {
      Log::misc() << "Evaluate Jacobian for Branch..." << std::endl;
      Log::misc() << "Jacobian evaluation is not implemented!" << std::endl;

      return 0;
    }

    // Available template instantiations
    template class Branch<DependencyTracking::Variable, long int>;
    template class Branch<DependencyTracking::Variable, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
