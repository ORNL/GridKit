#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>

#include "SampledSignalSourceImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace SignalSource
    {
      template class SampledSignalSource<DependencyTracking::Variable, long int>;
      template class SampledSignalSource<DependencyTracking::Variable, size_t>;
    } // namespace SignalSource
  } // namespace PhasorDynamics
} // namespace GridKit
