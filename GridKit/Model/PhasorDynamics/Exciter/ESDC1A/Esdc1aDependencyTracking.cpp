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
      template <class ScalarT, typename IdxT>
      int Esdc1a<ScalarT, IdxT>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for Esdc1a..." << std::endl;
        Log::misc() << "Jacobian evaluation not implemented!" << std::endl;
        return 0;
      }

      template class Esdc1a<DependencyTracking::Variable, long int>;
      template class Esdc1a<DependencyTracking::Variable, size_t>;
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
