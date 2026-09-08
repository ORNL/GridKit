#include "InnerCurrentControlImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      template class InnerCurrentControl<double, long int>;
      template class InnerCurrentControl<double, size_t>;
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
