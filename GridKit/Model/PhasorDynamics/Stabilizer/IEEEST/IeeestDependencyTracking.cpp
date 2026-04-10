/**
 * @file IeeestDependencyTracking.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief DependencyTracking Jacobian stub and template instantiations for IEEEST Stabilizer.
 */

#include "IeeestImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Stabilizer
    {
      /**
       * @brief Jacobian evaluation not implemented yet
       *
       * @tparam ScalarT - Scalar data type
       * @tparam IdxT    - Index data type
       * @return int - error code, 0 = success
       */
      template <class ScalarT, typename IdxT>
      int Ieeest<ScalarT, IdxT>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for Ieeest..." << std::endl;
        Log::misc() << "Jacobian evaluation not implemented!" << std::endl;
        return 0;
      }

      // Available template instantiations
      template class Ieeest<DependencyTracking::Variable, long int>;
      template class Ieeest<DependencyTracking::Variable, size_t>;
    } // namespace Stabilizer
  } // namespace PhasorDynamics
} // namespace GridKit
