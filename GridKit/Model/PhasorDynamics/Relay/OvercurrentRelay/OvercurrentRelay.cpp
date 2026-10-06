/**
 * @file OvercurrentRelay.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Non-Enzyme instantiation for the overcurrent relay model.
 */

#include "OvercurrentRelayImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Relay
    {
      /**
       * @brief Report that a separate Jacobian is unavailable in the plain-real build.
       */
      template <typename scalar_type, typename index_type>
      int OvercurrentRelay<scalar_type, index_type>::evaluateJacobian()
      {
        Log::misc() << "OvercurrentRelay: Jacobian evaluation is not implemented\n";
        return 0;
      }

      template class OvercurrentRelay<double, long int>;
      template class OvercurrentRelay<double, size_t>;
    } // namespace Relay
  } // namespace PhasorDynamics
} // namespace GridKit
