/**
 * @file GasBusInfinite.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 */

#include "GasBusInfiniteImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Gas
    {
      // Available template instantiations
      template class GasBusInfinite<double, long int>;
      template class GasBusInfinite<double, size_t>;
    } // namespace Gas
  } // namespace PhasorDynamics
} // namespace GridKit
