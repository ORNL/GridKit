
#include "GenClassicalImpl.hpp"

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
    int GenClassical<ScalarT, IdxT>::evaluateJacobian()
    {
      //std::cout << "Evaluate Jacobian for GenClassical..." << std::endl;
      //std::cout << "Jacobian evaluation is not implemented!" << std::endl;

      return 0;
    }

    // Available template instantiations
    template class GenClassical<double, long int>;
    template class GenClassical<double, size_t>;
  } // namespace PhasorDynamics
} // namespace GridKit
