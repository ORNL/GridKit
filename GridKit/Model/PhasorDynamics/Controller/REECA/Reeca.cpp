#include "ReecaImpl.hpp"

namespace GridKit::PhasorDynamics::Controller
{
  template <class ScalarT, class IdxT>
  int Reeca<ScalarT, IdxT>::evaluateJacobian()
  {
    return 0;
  }

  template class Reeca<double, long int>;
  template class Reeca<double, size_t>;
} // namespace GridKit::PhasorDynamics::Controller
