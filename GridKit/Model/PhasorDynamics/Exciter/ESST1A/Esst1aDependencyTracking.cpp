#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>

#include "Esst1aImpl.hpp"

namespace GridKit::PhasorDynamics::Exciter
{
  template <class ScalarT, class IdxT>
  int Esst1a<ScalarT, IdxT>::evaluateJacobian()
  {
    return this->constructCsr();
  }

  template class Esst1a<DependencyTracking::Variable, long int>;
  template class Esst1a<DependencyTracking::Variable, size_t>;
} // namespace GridKit::PhasorDynamics::Exciter
