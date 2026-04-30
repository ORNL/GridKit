#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>

#include "LoadImpl.hpp"

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
    int Load<ScalarT, IdxT>::evaluateJacobian()
    {
      Log::misc() << "Evaluate Jacobian for Load..." << std::endl;
      Log::misc() << "Jacobian evaluation is not implemented!" << std::endl;

      return 0;
    }

    // Available template instantiations
    template class Load<DependencyTracking::Variable, long int>;
    template class Load<DependencyTracking::Variable, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
