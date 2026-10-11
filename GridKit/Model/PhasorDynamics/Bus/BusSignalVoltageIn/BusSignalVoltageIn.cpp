/**
 * @file BusSignalVoltageIn.cpp
 * @author Slaven Peles (peless@ornl.gov)
 */

#include "BusSignalVoltageInImpl.hpp"

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
    int BusSignalVoltageIn<scalar_type, index_type>::evaluateJacobian()
    {
      Log::misc() << "BusSignalVoltageIn: Jacobian evaluation is not implemented\n";

      return 0;
    }

    // Available template instantiations
    template class BusSignalVoltageIn<double, long int>;
    template class BusSignalVoltageIn<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
