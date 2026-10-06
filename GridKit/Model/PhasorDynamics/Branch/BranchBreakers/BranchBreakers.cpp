/**
 * @file BranchBreakers.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Non-Enzyme instantiation for the breaker-terminated branch model.
 */

#include "BranchBreakersImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Report that a separate Jacobian is unavailable in the plain-real build.
     */
    template <typename scalar_type, typename index_type>
    int BranchBreakers<scalar_type, index_type>::evaluateJacobian()
    {
      Log::misc() << "BranchBreakers: Jacobian evaluation is not implemented\n";
      return 0;
    }

    template class BranchBreakers<double, long int>;
    template class BranchBreakers<double, size_t>;
  } // namespace PhasorDynamics
} // namespace GridKit
