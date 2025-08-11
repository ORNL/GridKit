#pragma once

#include <cstddef>
#include <type_traits>
#include <variant>
#include <vector>

#include "COO_Matrix.hpp"

namespace GridKit
{
  namespace LinearAlgebra
  {
    /**
     * @brief A helper mixin class which represents something (such as a matrix)
     * which may or may not be sorted.
     */
    class MaybeSorted
    {
    protected:
      /// @brief Flag which tracks whether or not the object is sorted
      bool sorted_;

      /**
       * @brief Cannot be constructed, except as a base class
       * @param isSorted Whether or not it should start sorted.
       */
      MaybeSorted(bool isSorted) : sorted_(isSorted)
      {
      }

    public:
      /**
       * @brief Whether or not it's currently sorted.
       */
      bool isSorted() const noexcept
      {
        return sorted_;
      }
    };

    /**
     * @brief A helper mixing class which (as opposed to `MaybeSorted`) represented something
     * (such as a matrix) in a constant state of being sorted.
     */
    class AlwaysSorted
    {
    protected:
      /**
       * @brief Cannot be constructed, except as a base class
       * @param isSorted Ignored - only for signature parity with `MaybeSorted`
       */
      AlwaysSorted([[maybe_unused]] bool isSorted)
      {
      }

    public:
      /**
       * @brief Whether or not it's currently sorted.
       * Always returns true, because it's always sorted.
       */
      constexpr bool isSorted() const noexcept
      {
        return true;
      }
    };

    /**
     * @brief Helper type - derive from this depending on whether your type is always sorted or not
     * @tparam KEEP_SORTED Whether your type is always sorted
     */
    template <bool KEEP_SORTED>
    using Sorted = std::conditional<KEEP_SORTED, AlwaysSorted, MaybeSorted>::type;

    template <class ScalarT, typename IdxT, bool KEEP_SORTED = false>
    class CSRMatrix : Sorted<KEEP_SORTED>
    {
    private:
      /**
       * @brief The array of (potentially) nonzero values.
       */
      std::vector<ScalarT> values_;
      /**
       * @brief The array of indices in the `::values_` array
       * that each row starts at. All values on a row exist purely
       * between that index (inclusive) and the next row's index (exclusive).
       * @note Some indices may be beyond the end of the `::values_` array -
       * this indicates that and all greater rows are all zero.
       */
      std::vector<IdxT>    row_indices_;
      /**
       * @brief The array of column indices for each value
       * in the `::values_` array. Each index must be strictly smaller than
       * `::num_cols_`.
       */
      std::vector<IdxT>    col_indices_;
      /**
       * @brief The number of columns of the matrix. The maximum value
       * (exclusive) for elements of `::col_indices_`.
       */
      IdxT                 num_cols_;

      CSRMatrix(std::vector<ScalarT> values, std::vector<IdxT> row_indices, std::vector<IdxT> col_indices, IdxT num_cols, bool sorted);

    public:
      CSRMatrix(size_t num_rows, IdxT num_cols);
      CSRMatrix(size_t num_rows, IdxT num_cols, size_t num_nonzero);

      static CSRMatrix fromCOO(COO_Matrix<ScalarT, IdxT>& coo);

      /**
       * @see `values_`
       */
      const std::vector<ScalarT>& values() const noexcept
      {
        return values_;
      }

      /**
       * @see `row_indices_`
       */
      const std::vector<IdxT>& rowIndices() const noexcept
      {
        return row_indices_;
      }

      /**
       * @see `col_indices_`
       */
      const std::vector<IdxT>& colIndices() const noexcept
      {
        return col_indices_;
      }

      /**
       * @brief The number of rows in the matrix. Expected to be
       * the size of `::row_indices_`.
       */
      size_t numRows() const noexcept
      {
        return row_indices_.size();
      }

      /**
       * @brief The number of columns in the matrix.
       */
      size_t numCols() const noexcept
      {
        return num_cols_;
      }

      /**
       * @brief Whether or not the column indices are currently sorted in each row.
       * If `KEEP_SORTED == true`, then this will always return `true`.
       * Only here for documentation purposes - this is really an inherited member.
       * @see sort()
       */
      bool isSorted() const noexcept
      {
        return Sorted<KEEP_SORTED>::isSorted();
      }

      void sort();
    };

    /**
     * @brief Construct a `CSRMatrix` from its parts.
     * @param sorted - Whether or not the `col_indices` array is sorted within each row.
     */
    template <class ScalarT, typename IdxT, bool KEEP_SORTED>
    CSRMatrix<ScalarT, IdxT, KEEP_SORTED>::CSRMatrix(std::vector<ScalarT> values, std::vector<IdxT> row_indices, std::vector<IdxT> col_indices, IdxT num_cols, bool sorted)
      : Sorted<KEEP_SORTED>(sorted), values_(values), row_indices_(row_indices), col_indices_(col_indices), num_cols_(num_cols)
    {
    }

    /**
     * @brief Create a new (empty) `CSRMatrix` of a specific size.
     */
    template <class ScalarT, typename IdxT, bool KEEP_SORTED>
    CSRMatrix<ScalarT, IdxT, KEEP_SORTED>::CSRMatrix(size_t num_rows, IdxT num_cols) : CSRMatrix(num_rows, num_cols, 0)
    {
    }

    /**
     * @brief Create a new (empty) `CSRMatrix` of a specific size,
     * and allocate enough space for a given amount of nonzero elements.
     */
    template <class ScalarT, typename IdxT, bool KEEP_SORTED>
    CSRMatrix<ScalarT, IdxT, KEEP_SORTED>::CSRMatrix(size_t num_rows, IdxT num_cols, size_t num_nonzero) : Sorted<KEEP_SORTED>(true), row_indices_(num_rows, 0), num_cols_(num_cols)
    {
      values_.reserve(num_nonzero);
      col_indices_.reserve(num_nonzero);
    }

    /**
     * @brief Construct a new `CSRMatrix` from the given `COO_Matrix`
     */
    template <class ScalarT, typename IdxT, bool KEEP_SORTED>
    CSRMatrix<ScalarT, IdxT, KEEP_SORTED> CSRMatrix<ScalarT, IdxT, KEEP_SORTED>::fromCOO(COO_Matrix<ScalarT, IdxT>& coo)
    {
      auto [row_indices, col_indices, values] = coo.setDataToCSR();
      auto [_num_rows, num_cols]              = coo.getDimensions();

      return CSRMatrix(values, row_indices, col_indices, num_cols, true);
    }
  } // namespace LinearAlgebra
} // namespace GridKit