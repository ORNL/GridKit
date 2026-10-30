/**
 * @file GasBranchDependencyTracking.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 */

#include "GasBranchImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Gas
    {
      /**
       * @brief Evaluate DependencyTracking::Variable Jacobian.
       *
       * @note Currently only used for testing.
       *
       * @return int - error code, 0 = success
       */
      template <typename scalar_type, typename index_type>
      int GasBranch<scalar_type, index_type>::evaluateJacobian()
      {
        this->constructCsr();

        return 0;
      }

      // Available template instantiations
      template class GasBranch<DependencyTracking::Variable, long int>;
      template class GasBranch<DependencyTracking::Variable, size_t>;
    } // namespace Gas
  } // namespace PhasorDynamics
} // namespace GridKit
