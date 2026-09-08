#include "OuterVoltageControlImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      template class OuterVoltageControl<double, long int>;
      template class OuterVoltageControl<double, size_t>;
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
