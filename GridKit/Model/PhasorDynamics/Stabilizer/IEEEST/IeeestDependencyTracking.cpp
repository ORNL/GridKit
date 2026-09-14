/**
 * @file IeeestDependencyTracking.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief DependencyTracking Jacobian stub and template instantiations for IEEEST Stabilizer.
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
      template <typename scalar_type, typename index_type, size_t order>
      int Ieeest<scalar_type, index_type, order>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for Ieeest..." << std::endl;
        Log::misc() << "Jacobian evaluation not implemented!" << std::endl;
        return 0;
      }

      // Available template instantiations
      template class Ieeest<DependencyTracking::Variable, long int, 0>;
      template class Ieeest<DependencyTracking::Variable, long int, 1>;
      template class Ieeest<DependencyTracking::Variable, long int, 2>;
      template class Ieeest<DependencyTracking::Variable, long int, 3>;
      template class Ieeest<DependencyTracking::Variable, long int, 4>;
      template class Ieeest<DependencyTracking::Variable, size_t, 0>;
      template class Ieeest<DependencyTracking::Variable, size_t, 1>;
      template class Ieeest<DependencyTracking::Variable, size_t, 2>;
      template class Ieeest<DependencyTracking::Variable, size_t, 3>;
      template class Ieeest<DependencyTracking::Variable, size_t, 4>;
    } // namespace Stabilizer
  } // namespace PhasorDynamics
} // namespace GridKit
