#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>

#include "Ieeeg1Impl.hpp"

namespace GridKit::PhasorDynamics::Governor
{
  template <class ScalarT, class IdxT>
  int Ieeeg1<ScalarT, IdxT>::evaluateJacobian()
  {
    return this->constructCsr();
  }

  template class Ieeeg1<DependencyTracking::Variable, long int>;
  template class Ieeeg1<DependencyTracking::Variable, size_t>;
} // namespace GridKit::PhasorDynamics::Governor
