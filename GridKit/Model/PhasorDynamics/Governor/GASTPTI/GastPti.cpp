/**
 * @file GastPti.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Non-Enzyme instantiation for the GASTPTI governor model.
 */

#include "GastPtiImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Governor
    {
      /**
       * @brief Report that a separate Jacobian is unavailable in the plain-real build.
       */
      template <typename scalar_type, typename index_type>
      int GastPti<scalar_type, index_type>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for GastPti..." << std::endl;
        Log::misc() << "Jacobian evaluation is not implemented!" << std::endl;
        return 0;
      }

      template class GastPti<double, long int>;
      template class GastPti<double, size_t>;
    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit
