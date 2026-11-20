/**
 * @file VoltageSourceDependencyTracking.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief EMT voltage source instantiation with dependency-tracking derivatives.
 */

#include <GridKit/Model/EMT/ComponentDependencyTracking.hpp>

#include "VoltageSourceImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    int VoltageSource<scalar_type, index_type>::evaluateResidual()
    {
      return ComponentT::evaluateResidual(*this);
    }

    template <typename scalar_type, typename index_type>
    int VoltageSource<scalar_type, index_type>::tagDifferentiable()
    {
      return ComponentT::tagDifferentiable(*this);
    }

    // Available template instantiations
    template class VoltageSource<double, size_t>;
  } // namespace EMT
} // namespace GridKit
