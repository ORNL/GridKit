/**
 * @file CoSimDependencyTracking.cpp
 * @author Philip Fackler (facklerpw@ornl.gov)
 *
 * @brief Deinition of CoSim connector interface.
 *
 */

#include "CoSimImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Connector
    {
      // Available template instantiations
      template class CoSim<DependencyTracking::Variable, long int>;
      template class CoSim<DependencyTracking::Variable, size_t>;

    } // namespace Connector
  } // namespace PhasorDynamics
} // namespace GridKit
