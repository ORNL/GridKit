
#include "BusImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Jacobian evaluation not implemented
     *
     * @tparam ScalarT - scalar data type
     * @tparam IdxT    - matrix index data type
     * @return int - error code, 0 = success
     */
    template <class ScalarT, typename IdxT>
    int Bus<ScalarT, IdxT>::evaluateJacobian()
    {
      std::cout << "Evaluate Jacobian for Bus..." << std::endl;
      std::cout << "Jacobian evaluation is not implemented!" << std::endl;

      return 0;
    }

    // Available template instantiations
    template class Bus<DependencyTracking::Variable, long int>;
    template class Bus<DependencyTracking::Variable, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
