#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>

#include "ReecaImpl.hpp"

namespace GridKit::PhasorDynamics::Controller
{
  template <class ScalarT, class IdxT>
  int Reeca<ScalarT, IdxT>::evaluateJacobian()
  {
    return this->constructCsr();
  }

  template class Reeca<DependencyTracking::Variable, long int>;
  template class Reeca<DependencyTracking::Variable, size_t>;
} // namespace GridKit::PhasorDynamics::Controller
