/**
 * @file GasBranch.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 */

#include "GasBranchImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Gas
    {
      /**
       * @brief Jacobian evaluation not implemented
       *
       * @return int - error code, 0 = success
       */
      template <typename scalar_type, typename index_type>
      int GasBranch<scalar_type, index_type>::evaluateJacobian()
      {
        Log::misc() << "GasBranch: Jacobian evaluation is not implemented\n";

        return 0;
      }

      // Available template instantiations
      template class GasBranch<double, long int>;
      template class GasBranch<double, size_t>;
    } // namespace Gas
  } // namespace PhasorDynamics
} // namespace GridKit
