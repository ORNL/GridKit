/**
 * @file BusSignalVoltageInDependencyTracking.cpp
 * @author Slaven Peles (peless@ornl.gov)
 */

#include "BusSignalVoltageInImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    // Available template instantiations
    template class BusSignalVoltageIn<DependencyTracking::Variable, long int>;
    template class BusSignalVoltageIn<DependencyTracking::Variable, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
