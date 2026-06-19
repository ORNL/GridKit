
#include "ConstantSignalSourceImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    // Available template instantiations
    template class ConstantSignalSource<DependencyTracking::Variable, long int>;
    template class ConstantSignalSource<DependencyTracking::Variable, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
