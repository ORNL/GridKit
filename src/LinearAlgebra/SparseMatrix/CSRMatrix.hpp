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
      MaybeSorted(bool isSorted)
        : sorted_(isSorted)
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
       * @note Some of the values may be zero, due to calculations working out that way.
       * This array stores exactly the values which *may* be nonzero.
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

      CSRMatrix(std::vector<ScalarT>&& values, std::vector<IdxT>&& row_indices, std::vector<IdxT>&& col_indices, IdxT num_cols, bool sorted);

    public:
      CSRMatrix(size_t num_rows, IdxT num_cols);
      CSRMatrix(size_t num_rows, IdxT num_cols, size_t num_nonzero);

      static CSRMatrix fromCOO(COO_Matrix<ScalarT, IdxT>& coo);
      static CSRMatrix fromCOO(COO_Matrix<ScalarT, IdxT>&& coo);

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
        return row_indices_.size() - 1;
      }

      /**
       * @brief The number of columns in the matrix.
       */
      size_t numCols() const noexcept
      {
        return num_cols_;
      }

      /**
       * @brief The number of (potentially) nonzero elements.
       */
      size_t numNonZero() const noexcept
      {
        return values_.size();
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

      ScalarT valueAt(size_t row, IdxT col) const;
    };

    /**
     * @brief Construct a `CSRMatrix` from its parts.
     * @param sorted - Whether or not the `col_indices` array is sorted within each row.
     */
    template <class ScalarT, typename IdxT, bool KEEP_SORTED>
    CSRMatrix<ScalarT, IdxT, KEEP_SORTED>::CSRMatrix(std::vector<ScalarT>&& values, std::vector<IdxT>&& row_indices, std::vector<IdxT>&& col_indices, IdxT num_cols, bool sorted)
      : Sorted<KEEP_SORTED>(sorted), values_(std::move(values)), row_indices_(std::move(row_indices)), col_indices_(std::move(col_indices)), num_cols_(num_cols)
    {
    }

    /**
     * @brief Create a new (empty) `CSRMatrix` of a specific size.
     */
    template <class ScalarT, typename IdxT, bool KEEP_SORTED>
    CSRMatrix<ScalarT, IdxT, KEEP_SORTED>::CSRMatrix(size_t num_rows, IdxT num_cols)
      : CSRMatrix(num_rows, num_cols, 0)
    {
    }

    /**
     * @brief Create a new (empty) `CSRMatrix` of a specific size,
     * and allocate enough space for a given amount of nonzero elements.
     */
    template <class ScalarT, typename IdxT, bool KEEP_SORTED>
    CSRMatrix<ScalarT, IdxT, KEEP_SORTED>::CSRMatrix(size_t num_rows, IdxT num_cols, size_t num_nonzero)
      : Sorted<KEEP_SORTED>(true), row_indices_(num_rows + 1, 0), num_cols_(num_cols)
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

      return CSRMatrix(std::move(values), std::move(row_indices), std::move(col_indices), std::move(num_cols), true);
    }

    /**
     * @brief Construct a new `CSRMatrix` from the given `COO_Matrix`. Faster with an r-value,
     * since we can repurpose `COO_Matrix::column_indices_` for our own `::col_indices_`.
     */
    template <class ScalarT, typename IdxT, bool KEEP_SORTED>
    CSRMatrix<ScalarT, IdxT, KEEP_SORTED> CSRMatrix<ScalarT, IdxT, KEEP_SORTED>::fromCOO(COO_Matrix<ScalarT, IdxT>&& coo)
    {
      auto row_indices           = coo.getCSRRowData();
      auto [_num_rows, num_cols] = coo.getDimensions();

      return CSRMatrix(std::move(coo.values_), std::move(row_indices), std::move(coo.column_indices_), std::move(num_cols), true);
    }

    /**
     * @brief Return the value of a particular element.
     * @param row The row of the element. This is bounds-checked in debug mode,
     *            but an out-of-bounds access otherwise is UB.
     * @param col The column of the elements. Bounds-checked always, since it needs to be searched for.
     * @return If the element exists in the matrix, returns its value. Otherwise return 0.
     * @note This access incurs a significant search cost, especially if the matrix is unsorted.
     */
    template <class ScalarT, typename IdxT, bool KEEP_SORTED>
    inline ScalarT CSRMatrix<ScalarT, IdxT, KEEP_SORTED>::valueAt(size_t row, IdxT col) const
    {
      // Bounds check for debug builds
      assert(row < numRows());

      // The indices in the values/col indices arrays that this row starts and ends at.
      // Our search for the element will be between these two spots
      IdxT row_start_idx = row_indices_[row];
      IdxT row_end_idx   = row_indices_[row + 1];

      auto row_start = col_indices_.begin() + row_start_idx;
      auto row_end   = col_indices_.begin() + row_end_idx;

      // The index in the values/col indices arrays we suspect the element may be found at.
      size_t found_idx;

      if (isSorted())
      {
        // Binary search can be used if the rows are sorted
        found_idx = std::lower_bound(row_start, row_end, col) - col_indices_.begin();
      }
      else
      {
        // If the rows aren't sorted, then we must do linear search
        found_idx = std::find(row_start, row_end, col) - col_indices_.begin();
      }

      // If the col we're trying to find wasn't found (either beyond the end of the row
      // or because we found a larger column index), the element isn't present in the
      // matrix, and we return 0.
      if (found_idx >= row_indices_[row + 1] || col_indices_[found_idx] != col)
      {
        return 0;
      }
      else
      {
        return values_[found_idx];
      }
    }
  } // namespace LinearAlgebra
} // namespace GridKit
