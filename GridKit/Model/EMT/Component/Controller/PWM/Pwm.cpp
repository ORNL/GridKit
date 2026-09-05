#include "PwmImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      template class Pwm<double, long int>;
      template class Pwm<double, size_t>;
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
