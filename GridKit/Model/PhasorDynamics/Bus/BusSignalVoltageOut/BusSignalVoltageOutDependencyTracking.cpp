/**
 * @file BusSignalVoltageOutDependencyTracking.cpp
 * @author Slaven Peles (peless@ornl.gov)
 */

#include "BusSignalVoltageOutImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Evaluate DependencyTracking::Variable Jacobian.
     *
     * DependencyTracking::Variable stores the Jacobian as dependency maps,
     * updated during calls to evaluateResidual(). Nothing to do here.
     *
     * @return int - error code
     */
    template <typename scalar_type, typename index_type>
    int BusSignalVoltageOut<scalar_type, index_type>::evaluateJacobian()
    {
      return 0;
    }

    // Available template instantiations
    template class BusSignalVoltageOut<DependencyTracking::Variable, long int>;
    template class BusSignalVoltageOut<DependencyTracking::Variable, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
