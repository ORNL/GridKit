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
    int Generator<scalar_type, index_type>::evaluateJacobian()
    {
      return ComponentT::evaluateJacobian(*this);
    }

    template <typename scalar_type, typename index_type>
    int Generator<scalar_type, index_type>::evaluateHessian(RealT sigma)
    {
      return ComponentT::evaluateHessian(*this, sigma);
    }

    // Available template instantiations
    template class Generator<double, size_t>;
  } // namespace OptimalPowerFlow
} // namespace GridKit
