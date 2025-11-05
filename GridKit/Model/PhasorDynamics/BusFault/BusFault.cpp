/**
 * @file BusFault.cpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Definition of a bus fault to ground model.
 *
 * The model uses Cartesian coordinates.
 */

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
      std::cout << "Evaluate Jacobian for BusFault..." << std::endl;
      std::cout << "Jacobian evaluation not implemented!" << std::endl;
      return 0;
    }

    // available template instantiations
    template class BusFault<double, long int>;
    template class BusFault<double, size_t>;
  } // namespace PhasorDynamics
} // namespace GridKit
