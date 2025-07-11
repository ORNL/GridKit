/**
 * @file BusFault.cpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Definition of a bus fault to ground model.
 *
 * The model uses Cartesian coordinates.
 */

#include "BusFaultImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    // available template instantiations
    template class BusFault<double, long int>;
    template class BusFault<double, size_t>;
  } // namespace PhasorDynamics
} // namespace GridKit
