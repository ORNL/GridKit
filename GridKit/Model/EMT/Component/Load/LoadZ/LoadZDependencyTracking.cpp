#include "LoadZImpl.hpp"

namespace GridKit::EMT
{
  template class LoadZ<DependencyTracking::Variable, long int>;
  template class LoadZ<DependencyTracking::Variable, size_t>;
} // namespace GridKit::EMT
