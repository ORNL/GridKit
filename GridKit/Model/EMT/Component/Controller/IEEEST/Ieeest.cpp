#include "IeeestJacobian.hpp"

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      template class Ieeest<double, long int>;
      template class Ieeest<double, size_t>;
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
