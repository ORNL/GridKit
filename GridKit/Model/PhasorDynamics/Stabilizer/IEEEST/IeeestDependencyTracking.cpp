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
       * @brief Evaluate DependencyTracking::Variable Jacobian.
       *
       * @note Currently only used for testing.
       *
       * DependencyTracking::Variable stores the Jacobian as dependency maps,
       * updated during calls to evaluateResidual().
       *
       * @return int - error code, 0 = success
       */
      template <typename scalar_type, typename index_type>
      int Ieeest<scalar_type, index_type>::evaluateJacobian()
      {
        Log::misc() << "Evaluate DependencyTracking Jacobian for Ieeest...\n";
        Log::misc() << "Jacobian evaluation is experimental!\n";

        this->constructCsr();

        return 0;
      }

      // Available template instantiations
      template class Ieeest<DependencyTracking::Variable, long int>;
      template class Ieeest<DependencyTracking::Variable, size_t>;

    } // namespace Stabilizer
  } // namespace PhasorDynamics
} // namespace GridKit
