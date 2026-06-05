/**
 * @file SexsPtiDependencyTracking.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Dependency-tracking instantiations for the SEXS-PTI exciter model.
 */

#include "SexsPtiImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      template <typename scalar_type, typename index_type>
      int SexsPti<scalar_type, index_type>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for SexsPti..." << std::endl;
        Log::misc() << "Jacobian evaluation not implemented!" << std::endl;
        return 0;
      }

      template class SexsPti<DependencyTracking::Variable, long int>;
      template class SexsPti<DependencyTracking::Variable, size_t>;

    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
