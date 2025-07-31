#include "Tgov1Impl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Governor
    {
      // Available template instantiations
      template class Tgov1<DependencyTracking::Variable, long int>;
      template class Tgov1<DependencyTracking::Variable, size_t>;
    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit
