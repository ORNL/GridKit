#include "GensalImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Jacobian evaluation not implemented yet
     *
     * @return int - error code, 0 = success
     */
    template <typename scalar_type, typename index_type>
    int Gensal<scalar_type, index_type>::evaluateJacobian()
    {
      Log::misc() << "Evaluate Jacobian for Gensal..." << std::endl;
      Log::misc() << "Jacobian evaluation not implemented!" << std::endl;
      return 0;
    }

    // Available template instantiations
    template class Gensal<DependencyTracking::Variable, long int>;
    template class Gensal<DependencyTracking::Variable, size_t>;
  } // namespace PhasorDynamics
} // namespace GridKit
