/**
 * @file Esdc1aEnzyme.cpp
 * @brief Optimal power flow ESDC1A exciter instantiation with Enzyme derivatives.
 */

#include <GridKit/Model/OptimalPowerFlow/ComponentEnzyme.hpp>

#include "Esdc1aImpl.hpp"

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    namespace Exciter
    {
      template <typename scalar_type, typename index_type>
      int Esdc1a<scalar_type, index_type>::evaluateGradient()
      {
        return ComponentT::evaluateGradient(*this);
      }

      template <typename scalar_type, typename index_type>
      int Esdc1a<scalar_type, index_type>::evaluateJacobian()
      {
        return ComponentT::evaluateJacobian(*this);
      }

      template <typename scalar_type, typename index_type>
      int Esdc1a<scalar_type, index_type>::evaluateHessian(RealT sigma)
      {
        return ComponentT::evaluateHessian(*this, sigma);
      }

      // Available template instantiations
      template class Esdc1a<double, size_t>;
    } // namespace Exciter
  } // namespace OptimalPowerFlow
} // namespace GridKit
