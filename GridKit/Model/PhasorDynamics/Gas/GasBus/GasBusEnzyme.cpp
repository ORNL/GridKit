/**
 * @file GasBusEnzyme.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 */

#include "GasBusImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Gas
    {
      /**
       * @brief The reset balance has no Jacobian entries; connected devices
       * add the entries of the bus row.
       *
       * @return int - error code
       */
      template <typename scalar_type, typename index_type>
      int GasBus<scalar_type, index_type>::evaluateJacobian()
      {
        return this->constructCoo();
      }

      // Available template instantiations
      template class GasBus<double, long int>;
      template class GasBus<double, size_t>;
    } // namespace Gas
  } // namespace PhasorDynamics
} // namespace GridKit
