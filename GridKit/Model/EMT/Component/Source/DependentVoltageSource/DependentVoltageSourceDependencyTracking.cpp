#include "DependentVoltageSourceImpl.hpp"

namespace GridKit::EMT
{
  template class DependentVoltageSource<DependencyTracking::Variable, long int>;
  template class DependentVoltageSource<DependencyTracking::Variable, size_t>;
} // namespace GridKit::EMT
