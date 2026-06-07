#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>

#include "DelayImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace SignalBlock
    {
      template class Delay<DependencyTracking::Variable, long int>;
      template class Delay<DependencyTracking::Variable, size_t>;
    } // namespace SignalBlock
  } // namespace PhasorDynamics
} // namespace GridKit
