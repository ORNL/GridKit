/**
 * @file GasCompressor.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 */

#include "GasCompressorImpl.hpp"

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
      int GasCompressor<scalar_type, index_type>::evaluateJacobian()
      {
        Log::misc() << "GasCompressor: Jacobian evaluation is not implemented\n";

        return 0;
      }

      // Available template instantiations
      template class GasCompressor<double, long int>;
      template class GasCompressor<double, size_t>;
    } // namespace Gas
  } // namespace PhasorDynamics
} // namespace GridKit
