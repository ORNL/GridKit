/**
 * @file Esdc1aDependencyTracking.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Dependency-tracking instantiations for the ESDC1A exciter model.
 */

#include "Esdc1aImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      /**
       * @brief Report that dependency tracking does not assemble a separate Jacobian.
       *
       * Dependency tracking exposes the Jacobian structure through the
       * residual rather than a separately assembled matrix.
       */
      template <typename scalar_type, typename index_type>
      int Esdc1a<scalar_type, index_type>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for Esdc1a..." << std::endl;
        Log::misc() << "Jacobian evaluation is not implemented!" << std::endl;
        return 0;
      }

      template class Esdc1a<DependencyTracking::Variable, long int>;
      template class Esdc1a<DependencyTracking::Variable, size_t>;
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
