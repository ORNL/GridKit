/**
 * @file RegcaDependencyTracking.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Dependency-tracking instantiations for the REGCA converter model.
 */

#include "RegcaImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Converter
    {
      // Dependency tracking recovers the sparsity pattern from the residual,
      // so this instantiation does not assemble a separate Jacobian.
      template <typename scalar_type, typename index_type>
      int Regca<scalar_type, index_type>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for Regca..." << std::endl;
        Log::misc() << "Jacobian evaluation is not implemented!" << std::endl;
        return 0;
      }

      template class Regca<DependencyTracking::Variable, long int>;
      template class Regca<DependencyTracking::Variable, size_t>;
    } // namespace Converter
  } // namespace PhasorDynamics
} // namespace GridKit
