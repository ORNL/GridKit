#include "DelaySmoothImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace SignalBlock
    {
      template <typename scalar_type, typename index_type>
      int DelaySmooth<scalar_type, index_type>::evaluateJacobian()
      {
        return 0;
      }

      template class DelaySmooth<DependencyTracking::Variable, long int>;
      template class DelaySmooth<DependencyTracking::Variable, size_t>;
    } // namespace SignalBlock
  } // namespace PhasorDynamics
} // namespace GridKit
