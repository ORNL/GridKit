/**
 * @file Bus.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief EMT bus instantiation.
 */

#include "BusImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    // Available template instantiations
    template class Bus<double, size_t>;
  } // namespace EMT
} // namespace GridKit
