
#include "FunctionSignalSourceImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    // Available template instantiations
    template class FunctionSignalSource<DependencyTracking::Variable, long int>;
    template class FunctionSignalSource<DependencyTracking::Variable, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
