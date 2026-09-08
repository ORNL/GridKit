#include "OuterVoltageControlImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      template class OuterVoltageControl<DependencyTracking::Variable, long int>;
      template class OuterVoltageControl<DependencyTracking::Variable, size_t>;
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
