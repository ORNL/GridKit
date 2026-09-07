#include "DcLinkImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      template class DcLink<DependencyTracking::Variable, long int>;
      template class DcLink<DependencyTracking::Variable, size_t>;
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
