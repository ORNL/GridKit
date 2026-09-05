#include "PwmImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      template class Pwm<DependencyTracking::Variable, long int>;
      template class Pwm<DependencyTracking::Variable, size_t>;
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
