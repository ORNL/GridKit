/**
 * @file GeneratorEnzyme.cpp
 * @brief Optimal power flow generator instantiation with Enzyme derivatives.
 */

#include <GridKit/AutomaticDifferentiation/Enzyme/D2LDx2.hpp>
#include <GridKit/AutomaticDifferentiation/Enzyme/DfDx.hpp>
#include <GridKit/AutomaticDifferentiation/Enzyme/DgDx.hpp>

#include "GeneratorImpl.hpp"

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    // Available template instantiations
    template class ComponentModel<Generator<double, size_t>, double, size_t>;
    template class Generator<double, size_t>;
  } // namespace OptimalPowerFlow
} // namespace GridKit
