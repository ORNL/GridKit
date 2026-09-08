#include "ModulationImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    template class Modulation<DependencyTracking::Variable, long int>;
    template class Modulation<DependencyTracking::Variable, size_t>;
  } // namespace EMT
} // namespace GridKit
