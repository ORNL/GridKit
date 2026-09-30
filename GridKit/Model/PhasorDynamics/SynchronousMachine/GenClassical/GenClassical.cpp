/**
 * @file GenClassical.cpp
 * @author Abdourahman Barry (abdourahman@vt.edu)
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Definition of a classical generator model.
 */

#include "GenClassicalImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Jacobian evaluation not implemented
     *
     * @return int - error code, 0 = success
     */
    template <typename scalar_type, typename index_type>
    int GenClassical<scalar_type, index_type>::evaluateJacobian()
    {
      Log::misc() << "GenClassical: Jacobian evaluation is not implemented\n";

      return 0;
    }

    // Available template instantiations
    template class GenClassical<double, long int>;
    template class GenClassical<double, size_t>;
  } // namespace PhasorDynamics
} // namespace GridKit
