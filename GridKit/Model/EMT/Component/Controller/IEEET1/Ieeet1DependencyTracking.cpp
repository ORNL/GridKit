#include "Ieeet1Jacobian.hpp"

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      template class Ieeet1<DependencyTracking::Variable, long int>;
      template class Ieeet1<DependencyTracking::Variable, size_t>;
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
