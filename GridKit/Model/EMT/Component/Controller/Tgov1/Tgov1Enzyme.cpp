/**
 * @file Tgov1Enzyme.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief EMT steam turbine governor instantiation with Enzyme derivatives.
 */

#include <GridKit/Model/EMT/ComponentEnzyme.hpp>

#include "Tgov1Impl.hpp"

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    int Tgov1<scalar_type, index_type>::evaluateResidual()
    {
      return ComponentT::evaluateResidual(*this);
    }

    template <typename scalar_type, typename index_type>
    int Tgov1<scalar_type, index_type>::evaluateJacobian()
    {
      return ComponentT::evaluateJacobian(*this);
    }

    template <typename scalar_type, typename index_type>
    int Tgov1<scalar_type, index_type>::tagDifferentiable()
    {
      return ComponentT::tagDifferentiable(*this);
    }

    // Available template instantiations
    template class Tgov1<double, size_t>;
  } // namespace EMT
} // namespace GridKit
