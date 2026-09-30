/**
 * @file Ieeet1.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @author Adam Birchfield (abirchfield@tamu.edu)
 *
 * @brief Definition of a IEEET1 Exciter.
 *
 */

#include "Ieeet1Impl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      /**
       * @brief Jacobian evaluation not implemented yet
       *
       * @return int - error code, 0 = success
       */
      template <typename scalar_type, typename index_type>
      int Ieeet1<scalar_type, index_type>::evaluateJacobian()
      {
        Log::misc() << "Ieeet1: Jacobian evaluation is not implemented\n";
        return 0;
      }

      // Available template instantiations
      template class Ieeet1<double, long int>;
      template class Ieeet1<double, size_t>;

    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
