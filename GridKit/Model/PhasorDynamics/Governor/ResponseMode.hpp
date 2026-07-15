/**
 * @file ResponseMode.hpp
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
      /// Governor response-limit behavior.
      enum class ResponseMode : std::size_t
      {
        Normal   = 0,
        Fixed    = 1,
        DownOnly = 2,
      };
    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit
