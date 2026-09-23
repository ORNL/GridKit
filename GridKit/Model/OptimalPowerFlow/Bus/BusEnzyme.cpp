/**
 * @file BusEnzyme.cpp
 * @brief Optimal power flow bus instantiation with Enzyme derivatives.
 */

#include <GridKit/AutomaticDifferentiation/Enzyme/D2LDx2.hpp>
#include <GridKit/AutomaticDifferentiation/Enzyme/DfDx.hpp>
#include <GridKit/AutomaticDifferentiation/Enzyme/DgDx.hpp>

#include "BusImpl.hpp"

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    // Available template instantiations
    template class ComponentModel<Bus<double, size_t>, double, size_t>;
    template class Bus<double, size_t>;
  } // namespace OptimalPowerFlow
} // namespace GridKit
