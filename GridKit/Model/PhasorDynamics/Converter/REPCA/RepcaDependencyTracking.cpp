/**
 * @file RepcaDependencyTracking.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Dependency-tracking instantiations for the REPCA plant-control model.
 */

#include "RepcaImpl.hpp"

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
      int Repca<scalar_type, index_type>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for Repca...\n";
        Log::misc() << "Jacobian evaluation is not implemented!\n";
        return 0;
      }

      template class Repca<DependencyTracking::Variable, long int>;
      template class Repca<DependencyTracking::Variable, size_t>;
    } // namespace Converter
  } // namespace PhasorDynamics
} // namespace GridKit
