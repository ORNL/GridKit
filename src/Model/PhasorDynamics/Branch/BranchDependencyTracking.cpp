/**
 * @file BranchDependencyTracking.cpp
 * @author Slaven Peles (peless@ornl.gov)
 *
 */

#include "BranchImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    // Available template instantiations
    template class Branch<DependencyTracking::Variable, long int>;
    template class Branch<DependencyTracking::Variable, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
