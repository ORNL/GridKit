/**
 * @file BusToSignalAdapterDependencyTracking.cpp
 * @author Philip Fackler (facklerpw@ornl.gov)
 *
 * @brief Deinition of BusToSignalAdapter connector interface.
 *
 */

#include "BusToSignalAdapterImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    // Available template instantiations
    template class BusToSignalAdapter<DependencyTracking::Variable, long int>;
    template class BusToSignalAdapter<DependencyTracking::Variable, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
