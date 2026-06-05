/**
 * @file GastPti.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Non-Enzyme instantiation for the GASTPTI governor model.
 */

#include "GastPtiImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Governor
    {
      template <class ScalarT, typename IdxT>
      int GastPti<ScalarT, IdxT>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for GastPti..." << std::endl;
        Log::misc() << "Jacobian evaluation is not implemented!" << std::endl;
        return 0;
      }

      template class GastPti<double, long int>;
      template class GastPti<double, size_t>;
    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit
