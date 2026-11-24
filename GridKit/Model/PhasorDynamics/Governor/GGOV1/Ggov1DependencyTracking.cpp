#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>

#include "Ggov1Impl.hpp"

namespace GridKit::PhasorDynamics::Governor
{
  template <class ScalarT, class IdxT>
  int Ggov1<ScalarT, IdxT>::evaluateJacobian()
  {
    return this->constructCsr();
  }

  template class Ggov1<DependencyTracking::Variable, long int>;
  template class Ggov1<DependencyTracking::Variable, size_t>;
} // namespace GridKit::PhasorDynamics::Governor
