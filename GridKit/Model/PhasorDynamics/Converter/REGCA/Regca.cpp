/**
 * @file Regca.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Non-Enzyme instantiation for the REGCA converter model.
 */

#include "RegcaImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Converter
    {
      /**
       * @brief Jacobian evaluation not implemented yet
       *
       * @return 0 on success.
       */
      template <typename scalar_type, typename index_type>
      int Regca<scalar_type, index_type>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for Regca..." << std::endl;
        Log::misc() << "Jacobian evaluation is not implemented!" << std::endl;
        return 0;
      }

      template class Regca<double, long int>;
      template class Regca<double, size_t>;
    } // namespace Converter
  } // namespace PhasorDynamics
} // namespace GridKit
