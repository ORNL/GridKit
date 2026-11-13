/**
 * @file RegcaEnzyme.cpp
 * @brief Optimal power flow REGCA converter instantiation with Enzyme derivatives.
 */

#include <GridKit/Model/OptimalPowerFlow/ComponentEnzyme.hpp>

#include "RegcaImpl.hpp"

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    namespace Converter
    {
      template <typename scalar_type, typename index_type>
      int Regca<scalar_type, index_type>::evaluateGradient()
      {
        return ComponentT::evaluateGradient(*this);
      }

      template <typename scalar_type, typename index_type>
      int Regca<scalar_type, index_type>::evaluateJacobian()
      {
        return ComponentT::evaluateJacobian(*this);
      }

      template <typename scalar_type, typename index_type>
      int Regca<scalar_type, index_type>::evaluateHessian(RealT sigma)
      {
        return ComponentT::evaluateHessian(*this, sigma);
      }

      // Available template instantiations
      template class Regca<double, size_t>;
    } // namespace Converter
  } // namespace OptimalPowerFlow
} // namespace GridKit
