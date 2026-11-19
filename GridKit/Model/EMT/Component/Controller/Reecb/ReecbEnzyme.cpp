/**
 * @file ReecbEnzyme.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief EMT renewable energy electrical control instantiation with Enzyme derivatives.
 */

#include <GridKit/Model/EMT/ComponentEnzyme.hpp>

#include "ReecbImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    int Reecb<scalar_type, index_type>::evaluateResidual()
    {
      return ComponentT::evaluateResidual(*this);
    }

    template <typename scalar_type, typename index_type>
    int Reecb<scalar_type, index_type>::evaluateJacobian()
    {
      return ComponentT::evaluateJacobian(*this);
    }

    template <typename scalar_type, typename index_type>
    int Reecb<scalar_type, index_type>::tagDifferentiable()
    {
      return ComponentT::tagDifferentiable(*this);
    }

    // Available template instantiations
    template class Reecb<double, size_t>;
  } // namespace EMT
} // namespace GridKit
