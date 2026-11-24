#include "Exst1PtiImpl.hpp"

namespace GridKit::PhasorDynamics::Exciter
{
  template <class ScalarT, class IdxT>
  int Exst1Pti<ScalarT, IdxT>::evaluateJacobian()
  {
    return 0;
  }

  template class Exst1Pti<double, long int>;
  template class Exst1Pti<double, size_t>;
} // namespace GridKit::PhasorDynamics::Exciter
