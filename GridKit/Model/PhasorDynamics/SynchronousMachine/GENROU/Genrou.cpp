/**
 * @file Genrou.cpp
 * @author Adam Birchfield (abirchfield@tamu.edu)
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Definition of a GENROU generator model.
 */

#include "GenrouImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Jacobian evaluation not implemented yet
     *
     * @return int - error code, 0 = success
     */
    template <typename scalar_type, typename index_type>
    int Genrou<scalar_type, index_type>::evaluateJacobian()
    {
      Log::misc() << "Genrou: Jacobian evaluation is not implemented\n";
      return 0;
    }

    // Available template instantiations
    template class Genrou<double, long int>;
    template class Genrou<double, size_t>;
  } // namespace PhasorDynamics
} // namespace GridKit
