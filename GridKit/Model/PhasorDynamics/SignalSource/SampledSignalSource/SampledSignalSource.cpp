#include "SampledSignalSourceImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace SignalSource
    {
      template class SampledSignalSource<double, long int>;
      template class SampledSignalSource<double, size_t>;
    } // namespace SignalSource
  } // namespace PhasorDynamics
} // namespace GridKit
