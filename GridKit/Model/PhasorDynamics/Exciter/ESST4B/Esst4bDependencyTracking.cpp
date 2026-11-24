#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>

#include "Esst4bImpl.hpp"

namespace GridKit::PhasorDynamics::Exciter
{
  template <class ScalarT, class IdxT>
  int Esst4b<ScalarT, IdxT>::evaluateJacobian()
  {
    return this->constructCsr();
  }

  template class Esst4b<DependencyTracking::Variable, long int>;
  template class Esst4b<DependencyTracking::Variable, size_t>;
} // namespace GridKit::PhasorDynamics::Exciter
