/**
 * @file Ieeet1Enzyme.cpp
 * @brief Optimal power flow IEEET1 exciter instantiation with Enzyme derivatives.
 */

#include <GridKit/Model/OptimalPowerFlow/ComponentEnzyme.hpp>

#include "Ieeet1Impl.hpp"

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    namespace Exciter
    {
      template <typename scalar_type, typename index_type>
      int Ieeet1<scalar_type, index_type>::evaluateGradient()
      {
        return ComponentT::evaluateGradient(*this);
      }

      template <typename scalar_type, typename index_type>
      int Ieeet1<scalar_type, index_type>::evaluateJacobian()
      {
        return ComponentT::evaluateJacobian(*this);
      }

      template <typename scalar_type, typename index_type>
      int Ieeet1<scalar_type, index_type>::evaluateHessian(RealT sigma)
      {
        return ComponentT::evaluateHessian(*this, sigma);
      }

      // Available template instantiations
      template class Ieeet1<double, size_t>;
    } // namespace Exciter
  } // namespace OptimalPowerFlow
} // namespace GridKit
