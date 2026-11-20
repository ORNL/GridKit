/**
 * @file BridgeEnzyme.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief EMT converter bridge instantiation with Enzyme derivatives.
 */

#include <GridKit/Model/EMT/ComponentEnzyme.hpp>

#include "BridgeImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    int Bridge<scalar_type, index_type>::evaluateResidual()
    {
      return evaluate([this](const auto& kernels)
                      { return this->ComponentT::evaluateResidual(kernels); });
    }

    template <typename scalar_type, typename index_type>
    int Bridge<scalar_type, index_type>::evaluateJacobian()
    {
      return evaluate([this](const auto& kernels)
                      { return this->ComponentT::evaluateJacobian(kernels); });
    }

    template <typename scalar_type, typename index_type>
    int Bridge<scalar_type, index_type>::tagDifferentiable()
    {
      return evaluate([this](const auto& kernels)
                      { return this->ComponentT::tagDifferentiable(kernels); });
    }

    // Available template instantiations
    template class Bridge<double, size_t>;
  } // namespace EMT
} // namespace GridKit
