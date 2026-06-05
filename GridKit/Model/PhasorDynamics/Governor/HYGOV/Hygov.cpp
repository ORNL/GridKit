/**
 * @file Hygov.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Non-Enzyme instantiation for the HYGOV governor model.
 */

#include "HygovImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Governor
    {
      template <class ScalarT, typename IdxT>
      int Hygov<ScalarT, IdxT>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for Hygov..." << std::endl;
        Log::misc() << "Jacobian evaluation is not implemented!" << std::endl;
        return 0;
      }

      template class Hygov<double, long int>;
      template class Hygov<double, size_t>;
    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit
