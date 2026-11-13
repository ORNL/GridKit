/**
 * @file HygovEnzyme.cpp
 * @brief Optimal power flow HYGOV governor instantiation with Enzyme derivatives.
 */

#include <GridKit/Model/OptimalPowerFlow/ComponentEnzyme.hpp>

#include "HygovImpl.hpp"

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    namespace Governor
    {
      template <typename scalar_type, typename index_type>
      int Hygov<scalar_type, index_type>::evaluateGradient()
      {
        return ComponentT::evaluateGradient(*this);
      }

      template <typename scalar_type, typename index_type>
      int Hygov<scalar_type, index_type>::evaluateJacobian()
      {
        return ComponentT::evaluateJacobian(*this);
      }

      // Available template instantiations
      template class Hygov<double, size_t>;
    } // namespace Governor
  } // namespace OptimalPowerFlow
} // namespace GridKit
