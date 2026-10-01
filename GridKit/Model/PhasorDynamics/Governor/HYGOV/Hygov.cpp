/**
 * @file Hygov.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Non-Enzyme instantiation for the HYGOV governor model.
 */

#include "HygovImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Governor
    {
      template <typename scalar_type, typename index_type>
      int Hygov<scalar_type, index_type>::evaluateJacobian()
      {
        Log::misc() << "Hygov: Jacobian evaluation is not implemented\n";
        return 0;
      }

      template class Hygov<double, long int>;
      template class Hygov<double, size_t>;
    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit
