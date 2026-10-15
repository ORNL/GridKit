/**
 * @file Load.cpp
 * @brief Optimal power flow load instantiation. Its kernels are constant, so
 * no Enzyme is needed.
 */

#include "LoadImpl.hpp"

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    // Available template instantiations
    template class ComponentModel<Load<double, size_t>, double, size_t>;
    template class Load<double, size_t>;
  } // namespace OptimalPowerFlow
} // namespace GridKit
