/**
 * @file RepcaEnzyme.cpp
 * @brief Optimal power flow REPCA controller instantiation with Enzyme derivatives.
 */

#include <GridKit/Model/OptimalPowerFlow/ComponentEnzyme.hpp>

#include "RepcaImpl.hpp"

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    namespace Controller
    {
      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::evaluateGradient()
      {
        return ComponentT::evaluateGradient(*this);
      }

      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::evaluateJacobian()
      {
        return ComponentT::evaluateJacobian(*this);
      }

      // Available template instantiations
      template class Repca<double, size_t>;
    } // namespace Controller
  } // namespace OptimalPowerFlow
} // namespace GridKit
