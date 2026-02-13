/**
 * @file Ieeet1.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @author Adam Birchfield (abirchfield@tamu.edu)
 *
 * @brief Definition of a IEEET1 Exciter.
 *
 */

#include "Ieeet1Impl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      /**
       * @brief Jacobian evaluation not implemented yet
       *
       * @tparam ScalarT - Scalar data type
       * @tparam IdxT    - Index data type
       * @return int - error code, 0 = success
       */
      template <class ScalarT, typename IdxT>
      int Ieeet1<ScalarT, IdxT>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for Ieeet1..." << std::endl;
        Log::misc() << "Jacobian evaluation not implemented!" << std::endl;
        return 0;
      }

      // Available template instantiations
      template class Ieeet1<double, long int>;
      template class Ieeet1<double, size_t>;

    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
