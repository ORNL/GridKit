#include "VoltageSourceImpl.hpp"

namespace GridKit::EMT
{
  template class VoltageSource<DependencyTracking::Variable, long int>;
  template class VoltageSource<DependencyTracking::Variable, size_t>;
} // namespace GridKit::EMT
