/**
 * @file PmuEnzyme.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief EMT phasor measurement unit instantiation with Enzyme derivatives.
 */

#include <GridKit/Model/EMT/ComponentEnzyme.hpp>

#include "PmuImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    int Pmu<scalar_type, index_type>::evaluateResidual()
    {
      return ComponentT::evaluateResidual(*this);
    }

    template <typename scalar_type, typename index_type>
    int Pmu<scalar_type, index_type>::evaluateOutputs()
    {
      return ComponentT::evaluateOutputs(*this);
    }

    template <typename scalar_type, typename index_type>
    int Pmu<scalar_type, index_type>::evaluateJacobian()
    {
      return ComponentT::evaluateJacobian(*this);
    }

    template <typename scalar_type, typename index_type>
    int Pmu<scalar_type, index_type>::tagDifferentiable()
    {
      return ComponentT::tagDifferentiable(*this);
    }

    // Available template instantiations
    template class Pmu<double, size_t>;
  } // namespace EMT
} // namespace GridKit
