/**
 * @file GastPtiEnzyme.cpp
 * @brief Optimal power flow GASTPTI governor instantiation with Enzyme derivatives.
 */

#include <GridKit/Model/OptimalPowerFlow/ComponentEnzyme.hpp>

#include "GastPtiImpl.hpp"

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    namespace Governor
    {
      template <typename scalar_type, typename index_type>
      int GastPti<scalar_type, index_type>::evaluateGradient()
      {
        return ComponentT::evaluateGradient(*this);
      }

      template <typename scalar_type, typename index_type>
      int GastPti<scalar_type, index_type>::evaluateJacobian()
      {
        return ComponentT::evaluateJacobian(*this);
      }

      // Available template instantiations
      template class GastPti<double, size_t>;
    } // namespace Governor
  } // namespace OptimalPowerFlow
} // namespace GridKit
