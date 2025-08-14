/**
 * @file Ieeet1DependencyTracking.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @author Adam Birchfield (abirchfield@tamu.edu)
 *
 * @brief Definition of a IEEET1 Exciter.
 *
 */

#include "Ieeet1Impl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {

      // Available template instantiations
      template class Ieeet1<DependencyTracking::Variable, long int>;
      template class Ieeet1<DependencyTracking::Variable, size_t>;

    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
