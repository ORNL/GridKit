/**
 * @file Tgov1DependencyTracking.cpp
 *
 */

#include "Tgov1Impl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Governor
    {
      /**
       * @brief Jacobian evaluation not implemented yet
       *
       * @tparam ScalarT - Scalar data type
       * @tparam IdxT    - Index data type
       * @return int - error code, 0 = success
       */
      template <class ScalarT, typename IdxT>
      int Tgov1<ScalarT, IdxT>::evaluateJacobian()
      {
        //std::cout << "Jacobian evaluation not implemented!" << std::endl;
        return 0;
      }

      // Available template instantiations
      template class Tgov1<DependencyTracking::Variable, long int>;
      template class Tgov1<DependencyTracking::Variable, size_t>;
    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit
