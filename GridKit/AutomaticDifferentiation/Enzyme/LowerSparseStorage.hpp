/**
 * @file LowerSparseStorage.hpp
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 *
 */

#pragma once

#include <vector>

namespace GridKit
{
  namespace Enzyme
  {
    namespace Sparse
    {
      /**
       * @brief Enzyme todense template
       *
       * @tparam T - return type
       */
      template <typename T>
      extern T __enzyme_todense(void*...) noexcept;

      /**
       * @brief Enzyme sparse storage in triplet format
       *
       * @tparam ScalarT - scalar data type
       */
      template <typename ScalarT>
      struct Triple
      {
        size_t  row;
        size_t  col;
        ScalarT val;
        Triple(Triple&&) = default;

        Triple(size_t row, size_t col, ScalarT val)
          : row(row),
            col(col),
            val(val)
        {
        }
      };

      /**
       * @brief Enzyme sparse accumulation for float
       *
       */
      [[maybe_unused]] __attribute__((enzyme_sparse_accumulate)) static void inner_storeflt(size_t row, size_t col, float val, std::vector<Triple<float>>& triplets)
      {
        triplets.emplace_back(row, col, val);
      }

      /**
       * @brief Enzyme sparse accumulation for double
       *
       */
      [[maybe_unused]] __attribute__((enzyme_sparse_accumulate)) static void inner_storedbl(size_t row, size_t col, double val, std::vector<Triple<double>>& triplets)
      {
        triplets.emplace_back(row, col, val);
      }

      /**
       * @brief Enzyme sparse store
       *
       * @tparam ScalarT - scalar data type
       */
      template <typename ScalarT>
      __attribute__((always_inline)) static void sparse_store(ScalarT val, size_t idx, size_t i, std::vector<Triple<ScalarT>>& triplets)
      {
        if (val == 0.0)
          return;
        idx /= sizeof(ScalarT);
        if constexpr (sizeof(ScalarT) == 4)
          inner_storeflt(idx, i, val, triplets);
        else
          inner_storedbl(idx, i, val, triplets);
      }

      /**
       * @brief Enzyme sparse load
       *
       * @tparam ScalarT - scalar data type
       */
      template <typename ScalarT>
      __attribute__((always_inline)) static ScalarT sparse_load(size_t, size_t, std::vector<Triple<ScalarT>>&)
      {
        return 0.0;
      }

      /**
       * @brief Enzyme identity store
       *
       * @tparam ScalarT - scalar data type
       */
      template <typename ScalarT>
      __attribute__((always_inline)) static void ident_store(ScalarT, size_t, size_t)
      {
        assert(0 && "should never store");
      }

      /**
       * @brief Enzyme identity load
       *
       * @tparam ScalarT - scalar data type
       */
      template <typename ScalarT>
      __attribute__((always_inline)) static ScalarT ident_load(size_t idx, size_t i)
      {
        idx /= sizeof(ScalarT);
        return (ScalarT) (idx == i);
      }
    } // namespace Sparse
  } // namespace Enzyme
} // namespace GridKit
