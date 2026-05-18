/**
 * @file ReecaDependencyTracking.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Dependency-tracking instantiations for the REECA converter-control model.
 */

#include "ReecaImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Converter
    {
      template <class ScalarT, typename IdxT>
      int Reeca<ScalarT, IdxT>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for Reeca..." << std::endl;
        Log::misc() << "Jacobian evaluation is not implemented!" << std::endl;
        return 0;
      }

      template class Reeca<DependencyTracking::Variable, long int>;
      template class Reeca<DependencyTracking::Variable, size_t>;
    } // namespace Converter
  } // namespace PhasorDynamics
} // namespace GridKit
