/**
 * @file InverseTimeRelay.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Non-Enzyme instantiation for the InverseTimeRelay model.
 */

#include "InverseTimeRelayImpl.hpp"

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
      int InverseTimeRelay<scalar_type, index_type>::evaluateJacobian()
      {
        Log::misc() << "InverseTimeRelay: Jacobian evaluation is not implemented\n";
        return 0;
      }

      template class InverseTimeRelay<double, long int>;
      template class InverseTimeRelay<double, size_t>;
    } // namespace Relay
  } // namespace PhasorDynamics
} // namespace GridKit
