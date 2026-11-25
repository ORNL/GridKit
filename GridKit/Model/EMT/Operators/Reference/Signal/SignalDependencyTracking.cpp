/**
 * @file SignalDependencyTracking.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief EMT signal instantiation with dependency-tracking derivatives.
 */

#include <GridKit/Model/EMT/ComponentDependencyTracking.hpp>

#include "SignalImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    int Signal<scalar_type, index_type>::tagDifferentiable()
    {
      return ComponentT::tagDifferentiable(*this);
    }

    // Available template instantiations
    template class Signal<double, size_t>;
  } // namespace EMT
} // namespace GridKit
