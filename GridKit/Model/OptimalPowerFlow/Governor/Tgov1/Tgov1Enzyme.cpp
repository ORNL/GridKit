/**
 * @file Tgov1Enzyme.cpp
 * @brief Optimal power flow TGOV1 governor instantiation with Enzyme derivatives.
 */

#include <GridKit/Model/OptimalPowerFlow/ComponentEnzyme.hpp>

#include "Tgov1Impl.hpp"

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    namespace Governor
    {
      template <typename scalar_type, typename index_type>
      int Tgov1<scalar_type, index_type>::evaluateGradient()
      {
        return ComponentT::evaluateGradient(*this);
      }

      template <typename scalar_type, typename index_type>
      int Tgov1<scalar_type, index_type>::evaluateJacobian()
      {
        return ComponentT::evaluateJacobian(*this);
      }

      // Available template instantiations
      template class Tgov1<double, size_t>;
    } // namespace Governor
  } // namespace OptimalPowerFlow
} // namespace GridKit
