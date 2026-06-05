/**
 * @file Ieeest.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Jacobian stub and template instantiations for IEEEST Stabilizer.
 */

#include "IeeestImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Stabilizer
    {
      /**
       * @brief Jacobian evaluation not implemented yet
       *
       * @return int - error code, 0 = success
       */
      template <typename scalar_type, typename index_type>
      int Ieeest<scalar_type, index_type>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for Ieeest..." << std::endl;
        Log::misc() << "Jacobian evaluation not implemented!" << std::endl;
        return 0;
      }

      // Available template instantiations
      template class Ieeest<double, long int>;
      template class Ieeest<double, size_t>;
    } // namespace Stabilizer
  } // namespace PhasorDynamics
} // namespace GridKit
