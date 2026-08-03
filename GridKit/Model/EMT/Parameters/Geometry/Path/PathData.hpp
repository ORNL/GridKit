/**
 * @file PathData.hpp
 *
 * @brief Static path data for finite-line EMT parameter models.
 *
 */

#pragma once

#include <optional>
#include <vector>

namespace GridKit
{
  namespace EMT
  {
    namespace Parameters
    {
      template <typename scalar_type>
      struct PathPoint
      {
        using ScalarT = scalar_type;

        ScalarT latitude{0.0};
        ScalarT longitude{0.0};
      };

      template <typename scalar_type, typename index_type>
      struct PathData
      {
        using ScalarT = scalar_type;
        using IdxT    = index_type;

        std::optional<ScalarT>          length;
        std::vector<PathPoint<ScalarT>> path;
      };
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
