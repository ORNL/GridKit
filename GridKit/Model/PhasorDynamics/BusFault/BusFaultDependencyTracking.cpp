#include "BusFaultImpl.hpp"

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
    int BusFault<scalar_type, index_type>::evaluateJacobian()
    {
      Log::misc() << "Evaluate Jacobian for BusFault..." << std::endl;
      Log::misc() << "Jacobian evaluation not implemented!" << std::endl;
      return 0;
    }

    // Additional template instantiations
    template class BusFault<DependencyTracking::Variable, long int>;
    template class BusFault<DependencyTracking::Variable, size_t>;
  } // namespace PhasorDynamics
} // namespace GridKit
