#include "SystemModelImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    // Available template instantiations
    // template class SystemModel<DependencyTracking::Variable, long int>;
    template class SystemModel<DependencyTracking::Variable, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
