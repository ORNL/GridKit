#include "GenrouImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    // Available template instantiations
    template class Genrou<DependencyTracking::Variable, long int>;
    template class Genrou<DependencyTracking::Variable, size_t>;
  } // namespace PhasorDynamics
} // namespace GridKit
