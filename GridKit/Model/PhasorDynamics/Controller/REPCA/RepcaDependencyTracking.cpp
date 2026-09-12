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
    namespace Controller
    {
      /**
       * @brief Evaluate DependencyTracking::Variable Jacobian.
       *
       * @note Currently only used for testing.
       *
       * DependencyTracking::Variable stores the Jacobian as dependency maps,
       * updated during calls to evaluateResidual().
       */
      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::evaluateJacobian()
      {
        Log::misc() << "Evaluate DependencyTracking Jacobian for Repca...\n";
        Log::misc() << "Jacobian evaluation is experimental!\n";

        this->constructCsr();

        return 0;
      }

      // Available template instantiations
      template class Repca<DependencyTracking::Variable, long int>;
      template class Repca<DependencyTracking::Variable, size_t>;

    } // namespace Controller
  } // namespace PhasorDynamics
} // namespace GridKit
