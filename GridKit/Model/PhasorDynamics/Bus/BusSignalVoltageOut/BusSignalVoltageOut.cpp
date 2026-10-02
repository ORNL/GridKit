/**
 * @file BusSignalVoltageOut.cpp
 * @author Slaven Peles (peless@ornl.gov)
 */

#include "BusSignalVoltageOutImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Jacobian evaluation not implemented
     *
     * @return int - error code
     */
    template <typename scalar_type, typename index_type>
    int BusSignalVoltageOut<scalar_type, index_type>::evaluateJacobian()
    {
      Log::misc() << "BusSignalVoltageOut: Jacobian evaluation is not implemented\n";

      return 0;
    }

    // Available template instantiations
    template class BusSignalVoltageOut<double, long int>;
    template class BusSignalVoltageOut<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
