#include "SexsPtiJacobian.hpp"

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      template class SexsPti<double, long int>;
      template class SexsPti<double, size_t>;
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
