/**
 * @file BusSignalVoltageIn.cpp
 * @author Slaven Peles (peless@ornl.gov)
 */

#include "BusSignalVoltageInImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    // Available template instantiations
    template class BusSignalVoltageIn<double, long int>;
    template class BusSignalVoltageIn<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
