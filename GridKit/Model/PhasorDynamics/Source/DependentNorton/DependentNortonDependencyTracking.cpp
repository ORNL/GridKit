/**
 * @file DependentNortonDependencyTracking.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Dependency-tracking instantiations for the dependent Norton source model.
 */

#include "DependentNortonImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Source
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
      int DependentNorton<scalar_type, index_type>::evaluateJacobian()
      {
        this->constructCsr();

        return 0;
      }

      // Available template instantiations
      template class DependentNorton<DependencyTracking::Variable, long int>;
      template class DependentNorton<DependencyTracking::Variable, size_t>;

    } // namespace Source
  } // namespace PhasorDynamics
} // namespace GridKit
