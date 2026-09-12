/**
 * @file SexsPtiDependencyTracking.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Dependency-tracking instantiations for the SEXS-PTI exciter model.
 */

#include "SexsPtiImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
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
      int SexsPti<scalar_type, index_type>::evaluateJacobian()
      {
        Log::misc() << "Evaluate DependencyTracking Jacobian for SexsPti...\n";
        Log::misc() << "Jacobian evaluation is experimental!\n";

        this->constructCsr();

        return 0;
      }

      // Available template instantiations
      template class SexsPti<DependencyTracking::Variable, long int>;
      template class SexsPti<DependencyTracking::Variable, size_t>;

    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
