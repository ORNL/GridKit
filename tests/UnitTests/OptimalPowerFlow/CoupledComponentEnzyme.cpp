#include <GridKit/Model/OptimalPowerFlow/ComponentEnzyme.hpp>

#include "CoupledComponent.hpp"

namespace GridKit
{
  namespace Testing
  {
    int CoupledComponent::evaluateGradient()
    {
      return Component::evaluateGradient(*this);
    }

    int CoupledComponent::evaluateEqualityJacobian()
    {
      return Component::evaluateEqualityJacobian(*this);
    }

    int CoupledComponent::evaluateInequalityJacobian()
    {
      return Component::evaluateInequalityJacobian(*this);
    }

    int CoupledComponent::evaluateHessian(RealT sigma, const RealT* lambda_c, const RealT* lambda_d)
    {
      return Component::evaluateHessian(*this, sigma, lambda_c, lambda_d);
    }
  } // namespace Testing
} // namespace GridKit
