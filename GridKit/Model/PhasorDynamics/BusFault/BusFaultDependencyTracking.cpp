#include "BusFaultImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Jacobian evaluation not implemented yet
     *
     * @tparam ScalarT - scalar data type
     * @tparam IdxT    - matrix index data type
     * @return int - error code, 0 = success
     */
    template <class ScalarT, typename IdxT>
    int BusFault<ScalarT, IdxT>::evaluateJacobian()
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
