/**
 * @file TransformerEnzyme.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief EMT transformer instantiation with Enzyme derivatives.
 */

#include <GridKit/Model/EMT/ComponentEnzyme.hpp>

#include "TransformerImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    int Transformer<scalar_type, index_type>::evaluateResidual()
    {
      return ComponentT::evaluateResidual(*this);
    }

    template <typename scalar_type, typename index_type>
    int Transformer<scalar_type, index_type>::evaluateJacobian()
    {
      return ComponentT::evaluateJacobian(*this);
    }

    template <typename scalar_type, typename index_type>
    int Transformer<scalar_type, index_type>::evaluateOutputs()
    {
      return ComponentT::evaluateOutputs(*this);
    }

    template <typename scalar_type, typename index_type>
    int Transformer<scalar_type, index_type>::tagDifferentiable()
    {
      return ComponentT::tagDifferentiable(*this);
    }

    // Available template instantiations
    template class Transformer<double, size_t>;
  } // namespace EMT
} // namespace GridKit
