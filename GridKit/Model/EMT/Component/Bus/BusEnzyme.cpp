/**
 * @file BusEnzyme.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 *
 */

#include "BusImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    // Available template instantiations
    template class Bus<double, long int>;
    template class Bus<double, size_t>;

  } // namespace EMT
} // namespace GridKit
