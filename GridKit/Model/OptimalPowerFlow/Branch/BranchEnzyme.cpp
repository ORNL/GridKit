/**
 * @file BranchEnzyme.cpp
 * @brief Optimal power flow branch instantiation with Enzyme derivatives.
 */

#include <GridKit/AutomaticDifferentiation/Enzyme/D2LDx2.hpp>
#include <GridKit/AutomaticDifferentiation/Enzyme/DfDx.hpp>
#include <GridKit/AutomaticDifferentiation/Enzyme/DgDx.hpp>

#include "BranchImpl.hpp"

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    // Available template instantiations
    template class ComponentModel<Branch<double, size_t>, double, size_t>;
    template class Branch<double, size_t>;
  } // namespace OptimalPowerFlow
} // namespace GridKit
