#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>

#include "Exst1PtiImpl.hpp"

namespace GridKit::PhasorDynamics::Exciter
{
  template <class ScalarT, class IdxT>
  int Exst1Pti<ScalarT, IdxT>::evaluateJacobian()
  {
    return this->constructCsr();
  }

  template class Exst1Pti<DependencyTracking::Variable, long int>;
  template class Exst1Pti<DependencyTracking::Variable, size_t>;
} // namespace GridKit::PhasorDynamics::Exciter
