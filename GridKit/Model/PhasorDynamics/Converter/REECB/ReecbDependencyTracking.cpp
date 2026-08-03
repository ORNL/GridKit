/**
 * @file ReecbDependencyTracking.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Dependency-tracking instantiations for the REECB electrical-control model.
 */

#include "ReecbImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Converter
    {
      /**
       * @brief Report that DependencyTracking exposes structure through the
       *        residual rather than a separately assembled Jacobian.
       */
      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for Reecb...\n";
        Log::misc() << "Jacobian evaluation is not implemented!\n";
        return 0;
      }

      template class Reecb<DependencyTracking::Variable, long int>;
      template class Reecb<DependencyTracking::Variable, size_t>;
    } // namespace Converter
  } // namespace PhasorDynamics
} // namespace GridKit
