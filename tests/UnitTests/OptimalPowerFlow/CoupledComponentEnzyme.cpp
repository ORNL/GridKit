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

    int CoupledComponent::evaluateJacobian()
    {
      return Component::evaluateJacobian(*this);
    }

    int CoupledComponent::evaluateHessian(RealT sigma)
    {
      return Component::evaluateHessian(*this, sigma);
    }
  } // namespace Testing
} // namespace GridKit
