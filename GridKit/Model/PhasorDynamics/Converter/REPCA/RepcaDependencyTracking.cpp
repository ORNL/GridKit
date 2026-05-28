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
      template <class ScalarT, typename IdxT>
      int Repca<ScalarT, IdxT>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for Repca..." << std::endl;
        Log::misc() << "Jacobian evaluation is not implemented!" << std::endl;
        return 0;
      }

      template class Repca<DependencyTracking::Variable, long int>;
      template class Repca<DependencyTracking::Variable, size_t>;
    } // namespace Converter
  } // namespace PhasorDynamics
} // namespace GridKit
