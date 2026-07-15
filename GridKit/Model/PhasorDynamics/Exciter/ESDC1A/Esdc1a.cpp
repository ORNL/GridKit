/**
 * @file Esdc1a.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Non-Enzyme instantiation for the ESDC1A exciter model.
 */

#include "Esdc1aImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      template <class ScalarT, typename IdxT>
      int Esdc1a<ScalarT, IdxT>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for Esdc1a..." << std::endl;
        Log::misc() << "Jacobian evaluation not implemented!" << std::endl;
        return 0;
      }

      template class Esdc1a<double, long int>;
      template class Esdc1a<double, size_t>;
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
