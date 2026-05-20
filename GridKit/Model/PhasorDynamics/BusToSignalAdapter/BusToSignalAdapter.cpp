/**
 * @file BusToSignalAdapter.cpp
 * @author Philip Fackler (facklerpw@ornl.gov)
 *
 * @brief Definition of BusToSignalAdapter component
 *
 */

#include "BusToSignalAdapterImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    // Available template instantiations
    template class BusToSignalAdapter<double, long int>;
    template class BusToSignalAdapter<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
