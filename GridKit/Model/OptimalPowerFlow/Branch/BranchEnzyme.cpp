/**
 * @file BranchEnzyme.cpp
 * @brief Optimal power flow branch instantiation with Enzyme derivatives.
 */

#include <GridKit/Model/OptimalPowerFlow/ComponentEnzyme.hpp>

#include "BranchImpl.hpp"

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    template <typename scalar_type, typename index_type>
    int Branch<scalar_type, index_type>::evaluateEqualityJacobian()
    {
      return ComponentT::evaluateEqualityJacobian(*this);
    }

    template <typename scalar_type, typename index_type>
    int Branch<scalar_type, index_type>::evaluateInequalityJacobian()
    {
      return ComponentT::evaluateInequalityJacobian(*this);
    }

    template <typename scalar_type, typename index_type>
    int Branch<scalar_type, index_type>::evaluateHessian(RealT sigma, const RealT* lambda_c, const RealT* lambda_d)
    {
      return ComponentT::evaluateHessian(*this, sigma, lambda_c, lambda_d);
    }

    // Available template instantiations
    template class Branch<double, size_t>;
  } // namespace OptimalPowerFlow
} // namespace GridKit
