/**
 * @file GasBus.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 */

#include "GasBusImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Gas
    {
      /**
       * @brief Jacobian evaluation not implemented
       *
       * @return int - error code
       */
      template <typename scalar_type, typename index_type>
      int GasBus<scalar_type, index_type>::evaluateJacobian()
      {
        Log::misc() << "GasBus: Jacobian evaluation is not implemented\n";

        return 0;
      }

      // Available template instantiations
      template class GasBus<double, long int>;
      template class GasBus<double, size_t>;
    } // namespace Gas
  } // namespace PhasorDynamics
} // namespace GridKit
