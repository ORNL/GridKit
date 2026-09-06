#include "GastPtiJacobian.hpp"

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      template class GastPti<double, long int>;
      template class GastPti<double, size_t>;
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
