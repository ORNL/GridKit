#include "InnerCurrentControlImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      template class InnerCurrentControl<DependencyTracking::Variable, long int>;
      template class InnerCurrentControl<DependencyTracking::Variable, size_t>;
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
