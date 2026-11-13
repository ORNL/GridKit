/**
 * @file ReecbEnzyme.cpp
 * @brief Optimal power flow REECB controller instantiation with Enzyme derivatives.
 */

#include <GridKit/Model/OptimalPowerFlow/ComponentEnzyme.hpp>

#include "ReecbImpl.hpp"

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    namespace Controller
    {
      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::evaluateGradient()
      {
        return ComponentT::evaluateGradient(*this);
      }

      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::evaluateJacobian()
      {
        return ComponentT::evaluateJacobian(*this);
      }

      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::evaluateHessian(RealT sigma)
      {
        return ComponentT::evaluateHessian(*this, sigma);
      }

      // Available template instantiations
      template class Reecb<double, size_t>;
    } // namespace Controller
  } // namespace OptimalPowerFlow
} // namespace GridKit
