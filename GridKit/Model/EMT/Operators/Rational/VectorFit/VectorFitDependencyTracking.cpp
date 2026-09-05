#include "VectorFit.hpp"

namespace GridKit
{
  namespace EMT
  {
    template class VectorFit<DependencyTracking::Variable, long int>;
    template class VectorFit<DependencyTracking::Variable, size_t>;
  } // namespace EMT
} // namespace GridKit
