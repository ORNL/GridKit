/**
 * @file HygovDependencyTracking.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Dependency-tracking instantiations for the HYGOV governor model.
 */

#include "HygovImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Governor
    {
      template <class ScalarT, typename IdxT>
      int Hygov<ScalarT, IdxT>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for Hygov..." << std::endl;
        Log::misc() << "Jacobian evaluation is not implemented!" << std::endl;
        return 0;
      }

      template class Hygov<DependencyTracking::Variable, long int>;
      template class Hygov<DependencyTracking::Variable, size_t>;
    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit
