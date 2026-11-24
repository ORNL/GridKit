#include "Ieeeg1Impl.hpp"

namespace GridKit::PhasorDynamics::Governor
{
  template <class ScalarT, class IdxT>
  int Ieeeg1<ScalarT, IdxT>::evaluateJacobian()
  {
    return 0;
  }

  template class Ieeeg1<double, long int>;
  template class Ieeeg1<double, size_t>;
} // namespace GridKit::PhasorDynamics::Governor
