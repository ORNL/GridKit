/**
 * @file ShuntEnzyme.cpp
 * @brief Optimal power flow shunt instantiation with Enzyme derivatives.
 */

#include <GridKit/AutomaticDifferentiation/Enzyme/D2LDx2.hpp>
#include <GridKit/AutomaticDifferentiation/Enzyme/DfDx.hpp>
#include <GridKit/AutomaticDifferentiation/Enzyme/DgDx.hpp>

#include "ShuntImpl.hpp"

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    // Available template instantiations
    template class ComponentModel<Shunt<double, size_t>, double, size_t>;
    template class Shunt<double, size_t>;
  } // namespace OptimalPowerFlow
} // namespace GridKit
