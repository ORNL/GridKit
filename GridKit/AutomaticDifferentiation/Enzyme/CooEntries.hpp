/**
 * @file CooEntries.hpp
 * @brief Coordinate entries filled by Enzyme sparse stores.
 */

#pragma once

#include <cstddef>
#include <vector>

namespace GridKit
{
  namespace Enzyme
  {
    namespace Sparse
    {
      /**
       * @brief Coordinate entries in global indices
       *
       * Filled by `mapped_store` and `lower_store`. Assembly sums duplicate
       * coordinates.
       */
      struct CooEntries
      {
        std::vector<size_t> rows;
        std::vector<size_t> cols;
        std::vector<double> values;

        void clear()
        {
          rows.clear();
          cols.clear();
          values.clear();
        }
      };
    } // namespace Sparse
  } // namespace Enzyme
} // namespace GridKit
