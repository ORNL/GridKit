/**
 * @file GastPtiDependencyTracking.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Dependency-tracking instantiations for the GASTPTI governor model.
 */

#include "GastPtiImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Governor
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
      int GastPti<scalar_type, index_type>::evaluateJacobian()
      {
        Log::misc() << "Evaluate DependencyTracking Jacobian for GastPti...\n";
        Log::misc() << "Jacobian evaluation is experimental!\n";

        this->constructCsr();

        return 0;
      }

      // Available template instantiations
      template class GastPti<DependencyTracking::Variable, long int>;
      template class GastPti<DependencyTracking::Variable, size_t>;

    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit
