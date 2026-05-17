/**
 * @file RegcaEnzyme.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Enzyme sparse Jacobian for the REGCA converter model.
 */

#include <GridKit/AutomaticDifferentiation/Enzyme/SparseJacobians.hpp>

#include "RegcaImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Converter
    {
      template <class ScalarT, typename IdxT>
      int Regca<ScalarT, IdxT>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for Regca..." << std::endl;
        Log::misc() << "Jacobian evaluation is experimental!" << std::endl;

        J_.zeroMatrix();
        return 0;
      }

      template class Regca<double, long int>;
      template class Regca<double, size_t>;
    } // namespace Converter
  } // namespace PhasorDynamics
} // namespace GridKit
