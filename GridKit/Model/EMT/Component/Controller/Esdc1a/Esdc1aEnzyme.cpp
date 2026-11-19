/**
 * @file Esdc1aEnzyme.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief EMT IEEE type DC1A excitation system instantiation with Enzyme derivatives.
 */

#include <GridKit/Model/EMT/ComponentEnzyme.hpp>

#include "Esdc1aImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    int Esdc1a<scalar_type, index_type>::evaluateResidual()
    {
      return ComponentT::evaluateResidual(*this);
    }

    template <typename scalar_type, typename index_type>
    int Esdc1a<scalar_type, index_type>::evaluateJacobian()
    {
      return ComponentT::evaluateJacobian(*this);
    }

    template <typename scalar_type, typename index_type>
    int Esdc1a<scalar_type, index_type>::tagDifferentiable()
    {
      return ComponentT::tagDifferentiable(*this);
    }

    // Available template instantiations
    template class Esdc1a<double, size_t>;
  } // namespace EMT
} // namespace GridKit
