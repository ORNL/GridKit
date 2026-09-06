#include "SexsPtiJacobian.hpp"

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      template class SexsPti<DependencyTracking::Variable, long int>;
      template class SexsPti<DependencyTracking::Variable, size_t>;
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
