/**
 * @file DependentNorton.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Real-scalar instantiations for the dependent Norton source model.
 */

#include "DependentNortonImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Source
    {
      /**
       * @brief Jacobian evaluation not implemented
       *
       * @return int - error code, 0 = success
       */
      template <typename scalar_type, typename index_type>
      int DependentNorton<scalar_type, index_type>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for DependentNorton..." << std::endl;
        Log::misc() << "Jacobian evaluation is not implemented!" << std::endl;

        return 0;
      }

      // Available template instantiations
      template class DependentNorton<double, long int>;
      template class DependentNorton<double, size_t>;

    } // namespace Source
  } // namespace PhasorDynamics
} // namespace GridKit
