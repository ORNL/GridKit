#include "GensalImpl.hpp"

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
    int Gensal<ScalarT, IdxT>::evaluateJacobian()
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
