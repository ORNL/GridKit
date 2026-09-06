#include "GastPtiJacobian.hpp"

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      template class GastPti<DependencyTracking::Variable, long int>;
      template class GastPti<DependencyTracking::Variable, size_t>;
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
