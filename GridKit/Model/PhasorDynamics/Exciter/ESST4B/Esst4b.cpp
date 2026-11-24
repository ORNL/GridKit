#include "Esst4bImpl.hpp"

namespace GridKit::PhasorDynamics::Exciter
{
  template <class ScalarT, class IdxT>
  int Esst4b<ScalarT, IdxT>::evaluateJacobian()
  {
    return 0;
  }

  template class Esst4b<double, long int>;
  template class Esst4b<double, size_t>;
} // namespace GridKit::PhasorDynamics::Exciter
