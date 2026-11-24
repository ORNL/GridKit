#include "Ggov1Impl.hpp"

namespace GridKit::PhasorDynamics::Governor
{
  template <class ScalarT, class IdxT>
  int Ggov1<ScalarT, IdxT>::evaluateJacobian()
  {
    return 0;
  }

  template class Ggov1<double, long int>;
  template class Ggov1<double, size_t>;
} // namespace GridKit::PhasorDynamics::Governor
