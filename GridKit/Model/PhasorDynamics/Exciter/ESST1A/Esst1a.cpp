#include "Esst1aImpl.hpp"

namespace GridKit::PhasorDynamics::Exciter
{
  template <class ScalarT, class IdxT>
  int Esst1a<ScalarT, IdxT>::evaluateJacobian()
  {
    return 0;
  }

  template class Esst1a<double, long int>;
  template class Esst1a<double, size_t>;
} // namespace GridKit::PhasorDynamics::Exciter
