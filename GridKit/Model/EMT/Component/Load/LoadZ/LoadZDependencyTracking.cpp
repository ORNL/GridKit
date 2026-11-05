/**
 * @file LoadZDependencyTracking.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief EMT impedance load instantiation with dependency-tracking derivatives.
 */

#include <GridKit/Model/EMT/ComponentDependencyTracking.hpp>

#include "LoadZImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    int LoadZ<scalar_type, index_type>::evaluateResidual()
    {
      return ComponentT::evaluateResidual(*this);
    }

    template <typename scalar_type, typename index_type>
    int LoadZ<scalar_type, index_type>::tagDifferentiable()
    {
      return ComponentT::tagDifferentiable(*this);
    }

    // Available template instantiations
    template class LoadZ<double, size_t>;
  } // namespace EMT
} // namespace GridKit
