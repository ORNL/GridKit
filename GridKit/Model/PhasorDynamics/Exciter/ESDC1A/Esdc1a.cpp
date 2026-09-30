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
      /**
       * @brief Report that a separate Jacobian is unavailable in the plain build.
       */
      template <typename scalar_type, typename index_type>
      int Esdc1a<scalar_type, index_type>::evaluateJacobian()
      {
        Log::misc() << "Esdc1a: Jacobian evaluation is not implemented\n";
        return 0;
      }

      template class Esdc1a<double, long int>;
      template class Esdc1a<double, size_t>;
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
