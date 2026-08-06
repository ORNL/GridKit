/**
 * @file Repca.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Non-Enzyme instantiation for the REPCA plant-control model.
 */

#include "RepcaImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Controller
    {
      /**
       * @brief Report that a separate Jacobian is unavailable in the plain-real build.
       */
      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for Repca...\n";
        Log::misc() << "Jacobian evaluation is not implemented!\n";
        return 0;
      }

      template class Repca<double, long int>;
      template class Repca<double, size_t>;
    } // namespace Controller
  } // namespace PhasorDynamics
} // namespace GridKit
