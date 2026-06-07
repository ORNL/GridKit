#include "DelayImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace SignalBlock
    {
      template class Delay<double, long int>;
      template class Delay<double, size_t>;
    } // namespace SignalBlock
  } // namespace PhasorDynamics
} // namespace GridKit
