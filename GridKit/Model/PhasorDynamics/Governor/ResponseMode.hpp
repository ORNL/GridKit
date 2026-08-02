/**
 * @file ResponseMode.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Governor response modes shared by PhasorDynamics governor models.
 */

#pragma once

#include <cstddef>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Governor
    {
      /// Governor response-limit behavior and serialized wire encoding.
      enum class ResponseMode : std::size_t
      {
        Normal   = 0, ///< \f$\mathrm{mode}=0\f$: use the configured response limits
        DownOnly = 1, ///< \f$\mathrm{mode}=1\f$: set the upper limit to the initial value
        Fixed    = 2, ///< \f$\mathrm{mode}=2\f$: set both limits to the initial value
      };
    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit
