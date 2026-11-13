/**
 * @file SexsPtiEnzyme.cpp
 * @brief Optimal power flow SEXS-PTI exciter instantiation with Enzyme derivatives.
 */

#include <GridKit/Model/OptimalPowerFlow/ComponentEnzyme.hpp>

#include "SexsPtiImpl.hpp"

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    namespace Exciter
    {
      template <typename scalar_type, typename index_type>
      int SexsPti<scalar_type, index_type>::evaluateGradient()
      {
        return ComponentT::evaluateGradient(*this);
      }

      template <typename scalar_type, typename index_type>
      int SexsPti<scalar_type, index_type>::evaluateJacobian()
      {
        return ComponentT::evaluateJacobian(*this);
      }

      // Available template instantiations
      template class SexsPti<double, size_t>;
    } // namespace Exciter
  } // namespace OptimalPowerFlow
} // namespace GridKit
