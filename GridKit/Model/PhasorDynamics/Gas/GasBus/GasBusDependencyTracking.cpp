/**
 * @file GasBusDependencyTracking.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 */

#include "GasBusImpl.hpp"

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
       * No-op, because connected devices track the bus-row dependencies.
       *
       * @return int - error code
       */
      template <typename scalar_type, typename index_type>
      int GasBus<scalar_type, index_type>::evaluateJacobian()
      {
        return 0;
      }

      // Available template instantiations
      template class GasBus<DependencyTracking::Variable, long int>;
      template class GasBus<DependencyTracking::Variable, size_t>;
    } // namespace Gas
  } // namespace PhasorDynamics
} // namespace GridKit
