#pragma once

#include <cstddef>
#include <numeric>
#include <optional>
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

    // Forward-declaration
    template <class, typename, bool, bool, bool>
    class CSRBuilder;

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

      // Explicit copying
      CSRMatrix(const CSRMatrix& other)            = default;
      CSRMatrix& operator=(const CSRMatrix& other) = default;

    public:
      CSRMatrix(CSRMatrix&& other) = default;
      CSRMatrix(size_t num_rows, IdxT num_cols);
      CSRMatrix(size_t num_rows, IdxT num_cols, size_t num_nonzero);

      CSRMatrix& operator=(CSRMatrix&& other) = default;

      static CSRMatrix fromCOO(COO_Matrix<ScalarT, IdxT>& coo);
      static CSRMatrix fromCOO(COO_Matrix<ScalarT, IdxT>&& coo);

      CSRMatrix clone() const;

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

      template <class, typename, bool, bool, bool>
      friend class CSRBuilder;
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
     * @brief Create a deep copy of this matrix. Since this operation is somewhat expensive,
     * it should be done explicitly.
     */
    template <class ScalarT, typename IdxT, bool KEEP_SORTED>
    CSRMatrix<ScalarT, IdxT, KEEP_SORTED> CSRMatrix<ScalarT, IdxT, KEEP_SORTED>::clone() const
    {
      return *this;
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

    /**
     * @brief A helpful builder utility to help build CSR matrices from their components. When built in debug mode,
     * enables several checks which ensure that you are building the CSR matrix correctly.
     * @tparam INCLUDE_DIAGONALS Whether or not we should ensure that empty diagonal elements are created, even if
     *                           not specified. This way if we axpy with a diagonal matrix later (a common operation
     *                           with Jacobian matrices) we can do it quickly.
     * @tparam KEEP_SORTED       Whether or not the sorting of rows should be enforced. If false, then elements in a row
     *                           can be inserted in arbitrary order.
     * @tparam USE_TEMPLATE      Whether or not a template matrix should be used. A template matrix should ideally be
     *                           one originally created using this utility, and whose internal buffers can be reused
     *                           for fast matrix construction.
     */
    template <class ScalarT, typename IdxT, bool INCLUDE_DIAGONALS = true, bool KEEP_SORTED = false, bool USE_TEMPLATE = true>
    class CSRBuilder
    {
      /**
       * @brief The matrix being built.
       */
      CSRMatrix<ScalarT, IdxT, KEEP_SORTED> mat_;
      /**
       * @brief The current row the builder is working on. Used to ensure we work on rows in ascending order,
       * and to insert extra diagonal elements if we skip rows.
       */
      size_t                                curr_row_       = 0;
      /**
       * @brief Keeps track of whether we have added a diagonal element in the row we're currently working on.
       * @todo  This can be removed from the class if `INCLUDE_DIAGONALS` is false. Maybe not necessary but an extra
       * little bit of performance/memory.
       */
      bool                                  added_diagonal_ = false;
      /**
       * @brief Keeps track of the next value to insert.
       * @todo  This can be removed from the class if `USE_TEMPLATE` is false.
       */
      size_t                                curr_val_       = 0;

      CSRBuilder(CSRMatrix<ScalarT, IdxT, KEEP_SORTED>&& mat) : mat_(std::move(mat))
      {
        if constexpr (!USE_TEMPLATE)
        {
          mat_.values_.clear();
          mat_.col_indices_.clear();
        }

        if constexpr (!KEEP_SORTED)
        {
          mat_.sorted_ = false;
        }
      }

      // Delete copy stuff - we're holding on to an expensive-to-copy matrix.
      CSRBuilder(const CSRBuilder& other)            = delete;
      CSRBuilder& operator=(const CSRBuilder& other) = delete;

      // (Weirdly) needed so that we can use constructors and such from instances of CSRBuilder with different
      // template values. See: fromEmpty and fromTemplate (which override USE_TEMPLATE). Not strictly necessary,
      // but this allows for a little bit of a better user experience.
      template <class, typename, bool, bool, bool>
      friend class CSRBuilder;

    public:
      CSRBuilder(CSRBuilder&& other)            = default;
      CSRBuilder& operator=(CSRBuilder&& other) = default;

      /**
       * @brief Start building a CSR matrix from an empty matrix of a specific size. Use this if you don't know what the final matrix
       * is going to look like, such as when you're constructing a Jacobian for the first time and can't use a previous Jacobian
       * as a template.
       * @param num_rows The number of rows of the matrix
       * @param num_cols The number of columns of the matrix
       * @return A new builder for the new matrix.
       */
      static CSRBuilder<ScalarT, IdxT, INCLUDE_DIAGONALS, KEEP_SORTED, false> fromEmpty(size_t num_rows, IdxT num_cols)
      {
        return CSRBuilder<ScalarT, IdxT, INCLUDE_DIAGONALS, KEEP_SORTED, false>(CSRMatrix<ScalarT, IdxT, KEEP_SORTED>(num_rows, num_cols));
      }

      /**
       * @brief Start building a CSR matrix from a template, which was prefereably constructed using a previous builder in the same way.
       * This will be faster than a builder constructed using `fromEmpty()`, since no re-allocation will need to be done, and in release mode
       * only values will be inserted (rows and columns will be ignored).
       * @param mat The matrix to use as a template.
       */
      static CSRBuilder<ScalarT, IdxT, INCLUDE_DIAGONALS, KEEP_SORTED, true> fromTemplate(CSRMatrix<ScalarT, IdxT, KEEP_SORTED>&& mat)
      {
        return CSRBuilder<ScalarT, IdxT, INCLUDE_DIAGONALS, KEEP_SORTED, true>(std::move(mat));
      }

      /**
       * @brief Finalize and return the constructed matrix.
       * @note  This irreparably finalizes the matrix being built - if an attempt is made to continue building the matrix after calling this,
       *        it will not work. In that way, we must move the builder to signal it is no longer usable afterwards.
       */
      operator CSRMatrix<ScalarT, IdxT, KEEP_SORTED>() &&
      {
        // We may have skipped some rows at the end - finish those entries up
        skipRows(mat_.numRows());

        // Setup the row index for the row beyond the last row
        nextRow(mat_.numRows());

        validate();
        return std::move(mat_);
      }

      /**
       * @brief Start constructing a row. Must be called before adding elements to the row.
       * @note Must be called in ascending order of rows. If `INCLUDE_DIAGONALS` is true, then this will add any skiped diagonal elements.
       * @return A reference to this builder, for builder chaining.
       */
      CSRBuilder& row(size_t new_row)
      {
        // Make sure this row is in-bounds
        assert(new_row < mat_.numRows());

        skipRows(new_row);

        if (new_row == 0)
        {
          // Ensure we haven't started adding values yet.
          assert(numNonZero() == 0);
        }
        else
        {
          // Otherwise ensure we're working in ascending order of rows
          assert(new_row > curr_row_);
          // Make sure all new values will be place correctly after the incoming row index
          assert(mat_.values_.size() > static_cast<size_t>(mat_.rowIndices()[new_row - 1]));
        }

        // Set new row
        nextRow(new_row);

        added_diagonal_ = false;

        return *this;
      }

      /**
       * @brief Insert a new element with a given value at the given column
       * @note If `KEEP_SORTED` is true, must be called in ascending order of column.
       * @return A reference to this builder, for builder chaining.
       */
      CSRBuilder& elem(IdxT column, ScalarT value)
      {
        if constexpr (INCLUDE_DIAGONALS)
        {
          // Insert diagonal element if we've gone past the diagonal
          if (!added_diagonal_ && static_cast<size_t>(column) > curr_row_)
          {
            if constexpr (KEEP_SORTED)
            {
              nextCol(static_cast<IdxT>(curr_row_));
              nextValue(0);

              added_diagonal_ = true;
            }
          }
          else
          {
            // We may be about to update the diagonal element
            added_diagonal_ = added_diagonal_ || (static_cast<size_t>(column) == curr_row_);
          }
        }

        if constexpr (KEEP_SORTED)
        {
          // Enforce ascending order of column
          assert(static_cast<size_t>(mat_.rowIndices()[curr_row_]) == numNonZero() || mat_.colIndices()[numNonZero() - 1] < column);
        }
        // Insert new column and value
        nextCol(column);
        nextValue(value);

        return *this;
      }

    private:
      /**
       * @brief Process the next column. Depending on the mode (`USE_TEMPLATE`), either inserts in column
       * or checks template and ensures the new column is as expected.
       */
      void nextCol(IdxT col)
      {
        if constexpr (USE_TEMPLATE)
        {
          assert(mat_.colIndices()[curr_val_] == col);
        }
        else
        {
          mat_.col_indices_.push_back(col);
        }
      }

      /**
       * @brief Process the next row. Depending on the mode (`USE_TEMPLATE`), either inserts
       * new row index, or checks template and ensures the row index is as expected.
       */
      void nextRow(size_t row)
      {
        curr_row_ = row;
        auto idx  = static_cast<IdxT>(numNonZero());

        if constexpr (USE_TEMPLATE)
        {
          assert(mat_.rowIndices()[curr_row_] == idx);
        }
        else
        {
          mat_.row_indices_[curr_row_] = idx;
        }
      }

      /**
       * @brief Inserts the next value. Depending on the mode (`USE_TEMPLATE`), either appends new value
       * to the end, or inserts values into where we expect it to go, and ensures we aren't overflowing the
       * template.
       * @param val
       */
      void nextValue(ScalarT val)
      {
        if constexpr (USE_TEMPLATE)
        {
          assert(curr_val_ < mat_.numNonZero());

          mat_.values_[curr_val_] = val;
          curr_val_++;
        }
        else
        {
          mat_.values_.push_back(val);
        }
      }

      /**
       * @brief The number of nonzero values we've inserted so far. Since we are re-using the
       * `CSRMatrix::values_` array when `USE_TEMPLATE` is true, `CSRMatrix::numNonZero()` will report
       * incorrectly.
       */
      size_t numNonZero() const noexcept
      {
        if constexpr (USE_TEMPLATE)
        {
          return curr_val_;
        }
        else
        {
          return mat_.numNonZero();
        }
      }

      /**
       * @brief A helper function for skipping rows. Sometimes a row gets skipped, such as when calling `row(1)` followed by `row(2)` or
       * by converting to a matrix before all rows have been entered in. This function will ensure the skipped rows properly have their
       * indices and diagonals added.
       */
      void skipRows(size_t until_row)
      {
        if constexpr (INCLUDE_DIAGONALS)
        {
          // Add the diagonal for the current row, if we missed it
          if (!added_diagonal_ && until_row > curr_row_)
          {
            elem(static_cast<IdxT>(curr_row_), 0);
          }
        }

        // Skipped rows
        for (size_t row = curr_row_ + 1; row < until_row; row++)
        {
          nextRow(row);

          // Add diagonals for skipped rows if needed
          if constexpr (INCLUDE_DIAGONALS)
          {
            elem(static_cast<IdxT>(row), 0);
          }
        }
      }

      /**
       * @brief Ensure that the matrix under construction is valid, and that, if using a template,
       * the template was used as expected.
       */
      void validate() const
      {
        // The first row should always be pointing to the first value - even if empty.
        assert(mat_.rowIndices()[0] == 0);

        for (size_t i = 1; i < mat_.numRows(); i++)
        {
          // Ascending order of rows
          assert(mat_.rowIndices()[i] >= mat_.rowIndices()[i - 1]);
          // In-bounds
          assert(static_cast<size_t>(mat_.rowIndices()[i]) <= numNonZero());

          bool check_diagonal = !INCLUDE_DIAGONALS;
          for (IdxT j = mat_.rowIndices()[i]; j < mat_.rowIndices()[i + 1]; j++)
          {
            if constexpr (KEEP_SORTED)
            {
              // Ascending order of columns
              if (j < mat_.rowIndices()[i + 1] - 1)
              {
                assert(mat_.colIndices()[j] < mat_.colIndices()[j + 1]);
              }
            }

            // In-bounds
            assert(static_cast<size_t>(mat_.colIndices()[j]) < mat_.numCols());

            // Check to make sure we have a diagonal element
            check_diagonal = check_diagonal || static_cast<size_t>(mat_.colIndices()[j]) == i;
          }

          assert(check_diagonal);
        }
        // Make sure last row index is just out of bounds
        assert(static_cast<size_t>(mat_.rowIndices().back()) == mat_.numNonZero());

        if constexpr (USE_TEMPLATE)
        {
          // Make sure we added as many elements as the template had
          assert(numNonZero() == mat_.colIndices().size());
          // Make sure we went through as many rows as the template had
          assert(curr_row_ == mat_.numRows());
        }

        // If we aren't ensuring the building is done in ascending order of column,
        // then there isn't anything checking that each column is present at most once
        // in each row. Check that here.
        if constexpr (!KEEP_SORTED)
        {
          std::vector<IdxT> col_indices = mat_.col_indices_;
          for (size_t i = 0; i < mat_.numRows(); i++)
          {
            std::sort(col_indices.begin() + mat_.row_indices_[i], col_indices.begin() + mat_.row_indices_[i + 1]);

            for (size_t j = mat_.row_indices_[i]; j < mat_.row_indices_[i + 1] - 1; j++)
            {
              assert(col_indices[j] != col_indices[j + 1]);
            }
          }
        }
      }
    };
  } // namespace LinearAlgebra
} // namespace GridKit
