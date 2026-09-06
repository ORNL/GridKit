#include "IeeestJacobian.hpp"

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      template class Ieeest<DependencyTracking::Variable, long int>;
      template class Ieeest<DependencyTracking::Variable, size_t>;
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
