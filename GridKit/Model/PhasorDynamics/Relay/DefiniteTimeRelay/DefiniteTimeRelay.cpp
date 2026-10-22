/**
 * @file DefiniteTimeRelay.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Non-Enzyme instantiation for the DefiniteTimeRelay model.
 */

#include "DefiniteTimeRelayImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Relay
    {
      /**
       * @brief Report that a separate Jacobian is unavailable in the plain build.
       */
      template <typename scalar_type, typename index_type>
      int DefiniteTimeRelay<scalar_type, index_type>::evaluateJacobian()
      {
        Log::misc() << "DefiniteTimeRelay: Jacobian evaluation is not implemented\n";
        return 0;
      }

      template class DefiniteTimeRelay<double, long int>;
      template class DefiniteTimeRelay<double, size_t>;
    } // namespace Relay
  } // namespace PhasorDynamics
} // namespace GridKit
