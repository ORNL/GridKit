/**
 * @file GeneratorEnzyme.cpp
 * @brief Optimal power flow generator instantiation with Enzyme derivatives.
 */

#include <GridKit/Model/OptimalPowerFlow/ComponentEnzyme.hpp>

#include "GeneratorImpl.hpp"

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    template <typename scalar_type, typename index_type>
    int Generator<scalar_type, index_type>::evaluateGradient()
    {
      return ComponentT::evaluateGradient(*this);
    }

    template <typename scalar_type, typename index_type>
    int Generator<scalar_type, index_type>::evaluateEqualityJacobian()
    {
      return ComponentT::evaluateEqualityJacobian(*this);
    }

    template <typename scalar_type, typename index_type>
    int Generator<scalar_type, index_type>::evaluateHessian(RealT sigma, const RealT* lambda_c, const RealT* lambda_d)
    {
      return ComponentT::evaluateHessian(*this, sigma, lambda_c, lambda_d);
    }

    // Available template instantiations
    template class Generator<double, size_t>;
  } // namespace OptimalPowerFlow
} // namespace GridKit
