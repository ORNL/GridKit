/**
 * @file TowerData.hpp
 *
 * @brief Static conductor tower data for overhead line parameter models.
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
      template <typename scalar_type, typename index_type>
      struct TowerData
      {
        using ScalarT = scalar_type;
        using IdxT    = index_type;

        IdxT                                K{0};
        ScalarT                             span{0.0};
        std::vector<ScalarT>                height;
        std::vector<ScalarT>                position;
        std::optional<std::vector<ScalarT>> tension;
      };
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
