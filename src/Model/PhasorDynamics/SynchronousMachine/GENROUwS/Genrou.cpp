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
    // Available template instantiations
    template class Genrou<double, long int>;
    template class Genrou<double, size_t>;
  } // namespace PhasorDynamics
} // namespace GridKit
