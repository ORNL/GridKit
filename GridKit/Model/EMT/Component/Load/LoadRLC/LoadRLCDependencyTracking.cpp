/**
 * @file LoadRLCDependencyTracking.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief EMT parallel RLC load instantiation with dependency-tracking derivatives.
 */

#include <GridKit/Model/EMT/ComponentDependencyTracking.hpp>

#include "LoadRLCImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    int LoadRLC<scalar_type, index_type>::evaluateOutputs()
    {
      return ComponentT::evaluateOutputs(*this);
    }

    template <typename scalar_type, typename index_type>
    int LoadRLC<scalar_type, index_type>::tagDifferentiable()
    {
      return ComponentT::tagDifferentiable(*this);
    }

    // Available template instantiations
    template class LoadRLC<double, size_t>;
  } // namespace EMT
} // namespace GridKit
