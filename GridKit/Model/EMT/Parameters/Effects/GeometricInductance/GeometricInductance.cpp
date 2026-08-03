/**
 * @file GeometricInductance.cpp
 *
 * @brief Sparse-jet Jacobian for the GeometricInductance model.
 *
 */

#include "GeometricInductance.hpp"

namespace GridKit
{
  namespace EMT
  {
    namespace Parameters
    {
      template <typename ScalarT, typename IdxT>
      int GeometricInductance<ScalarT, IdxT>::evaluateJacobian()
      {
        return this->template evaluateLocalAlgebraicJacobian<GeometricInductance<ScalarT, IdxT>>();
      }

      template class GeometricInductance<double, long int>;
      template class GeometricInductance<double, size_t>;
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
