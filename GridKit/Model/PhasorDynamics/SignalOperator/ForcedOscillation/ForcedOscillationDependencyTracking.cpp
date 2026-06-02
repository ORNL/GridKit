/**
 * @file ForcedOscillationDependencyTracking.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Dependency tracking instantiations for ForcedOscillation.
 */

#include "ForcedOscillationImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace SignalOperator
    {
      /**
       * @brief Jacobian evaluation not implemented yet
       *
       * @tparam ScalarT - Scalar data type
       * @tparam IdxT    - Index data type
       * @return int - error code, 0 = success
       */
      template <class ScalarT, typename IdxT>
      int ForcedOscillation<ScalarT, IdxT>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for ForcedOscillation..." << std::endl;
        Log::misc() << "Jacobian evaluation not implemented!" << std::endl;
        return 0;
      }

      // Available template instantiations
      template class ForcedOscillation<DependencyTracking::Variable, long int>;
      template class ForcedOscillation<DependencyTracking::Variable, size_t>;

    } // namespace SignalOperator
  } // namespace PhasorDynamics
} // namespace GridKit
