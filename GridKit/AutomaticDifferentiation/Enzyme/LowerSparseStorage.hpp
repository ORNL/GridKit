/**
 * @file LowerSparseStorage.hpp
 *
 * @details This file contains functions used by Enzyme to store sparse Jacobians.
 *
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 *
 */

#pragma once

#include <cassert>
#include <vector>

#include <GridKit/AutomaticDifferentiation/Enzyme/CooEntries.hpp>
#include <GridKit/Constants.hpp>

namespace GridKit
{
  namespace Enzyme
  {
    namespace Sparse
    {
      /**
       * @brief Enzyme todense template
       *
       * @details This is used by Enzyme's auto sparsity analysis. It internally maps dense storage
       * to sparse ones (where structural zeros are not kept).
       *
       * @tparam T - return type
       */
      template <typename T>
      extern T __enzyme_todense(void*...) noexcept;

      /**
       * @brief Enzyme sparse accumulation for float and size_t
       *
       * @note __attribute__((enzyme_sparse_accumulate)) does not support templates yet
       *
       * @param[in] row - row to be stored
       * @param[in] col - column to be stored
       * @param[in] val - value to be stored
       * @param[in] scaling - scaling factor for values
       * @param[in] res_indices - Global residual indices
       * @param[in] var_indices - Global variable indices
       * @param[in,out] rows - buffer where row will be stored
       * @param[in,out] cols - buffer where col will be stored
       * @param[in,out] vals - buffer where val will be stored
       * @param[in,out] nnz - number of nonzeros
       */
      [[maybe_unused]] __attribute__((enzyme_sparse_accumulate)) static void inner_store_float_size_t(
          size_t        row,
          size_t        col,
          float         val,
          float         scaling,
          const size_t* row_indices,
          const size_t* col_indices,
          size_t*       rows,
          size_t*       cols,
          float*        vals,
          size_t&       nnz)
      {
        const auto row_mapped = row_indices[static_cast<size_t>(row)];
        const auto col_mapped = col_indices[static_cast<size_t>(col)];
        if (col_mapped != INVALID_INDEX<size_t>)
        {
          rows[static_cast<size_t>(nnz)] = row_mapped;
          cols[static_cast<size_t>(nnz)] = col_mapped;
          vals[static_cast<size_t>(nnz)] = scaling * val;
          nnz++;
        }
      }

      /**
       * @brief Enzyme sparse accumulation for float and long int
       *
       * @note __attribute__((enzyme_sparse_accumulate)) does not support templates yet
       *
       * @param[in] row - row to be stored
       * @param[in] col - column to be stored
       * @param[in] val - value to be stored
       * @param[in] scaling - scaling factor for values
       * @param[in] res_indices - Global residual indices
       * @param[in] var_indices - Global variable indices
       * @param[in,out] rows - buffer where row will be stored
       * @param[in,out] cols - buffer where col will be stored
       * @param[in,out] vals - buffer where val will be stored
       * @param[in,out] nnz - number of nonzeros
       */
      [[maybe_unused]] __attribute__((enzyme_sparse_accumulate)) static void inner_store_float_long_int(
          long int        row,
          long int        col,
          float           val,
          float           scaling,
          const long int* row_indices,
          const long int* col_indices,
          long int*       rows,
          long int*       cols,
          float*          vals,
          long int&       nnz)
      {
        const auto row_mapped = row_indices[static_cast<size_t>(row)];
        const auto col_mapped = col_indices[static_cast<size_t>(col)];
        if (col_mapped != INVALID_INDEX<long int>)
        {
          rows[static_cast<size_t>(nnz)] = row_mapped;
          cols[static_cast<size_t>(nnz)] = col_mapped;
          vals[static_cast<size_t>(nnz)] = scaling * val;
          nnz++;
        }
      }

      /**
       * @brief Enzyme sparse accumulation for double and size_t
       *
       * @note __attribute__((enzyme_sparse_accumulate)) does not support templates yet
       *
       * @param[in] row - row to be stored
       * @param[in] col - column to be stored
       * @param[in] val - value to be stored
       * @param[in] scaling - scaling factor for values
       * @param[in] res_indices - Global residual indices
       * @param[in] var_indices - Global variable indices
       * @param[in,out] rows - buffer where row will be stored
       * @param[in,out] cols - buffer where col will be stored
       * @param[in,out] vals - buffer where val will be stored
       * @param[in,out] nnz - number of nonzeros
       */
      [[maybe_unused]] __attribute__((enzyme_sparse_accumulate)) static void inner_store_double_size_t(
          size_t        row,
          size_t        col,
          double        val,
          double        scaling,
          const size_t* row_indices,
          const size_t* col_indices,
          size_t*       rows,
          size_t*       cols,
          double*       vals,
          size_t&       nnz)
      {
        const auto row_mapped = row_indices[static_cast<size_t>(row)];
        const auto col_mapped = col_indices[static_cast<size_t>(col)];
        if (col_mapped != INVALID_INDEX<size_t>)
        {
          rows[static_cast<size_t>(nnz)] = row_mapped;
          cols[static_cast<size_t>(nnz)] = col_mapped;
          vals[static_cast<size_t>(nnz)] = scaling * val;
          nnz++;
        }
      }

      /**
       * @brief Enzyme sparse accumulation for double and long int
       *
       * @note __attribute__((enzyme_sparse_accumulate)) does not support templates yet
       *
       * @param[in] row - row to be stored
       * @param[in] col - column to be stored
       * @param[in] val - value to be stored
       * @param[in] scaling - scaling factor for values
       * @param[in] res_indices - Global residual indices
       * @param[in] var_indices - Global variable indices
       * @param[in,out] rows - buffer where row will be stored
       * @param[in,out] cols - buffer where col will be stored
       * @param[in,out] vals - buffer where val will be stored
       * @param[in,out] nnz - number of nonzeros
       */
      [[maybe_unused]] __attribute__((enzyme_sparse_accumulate)) static void inner_store_double_long_int(
          long int        row,
          long int        col,
          double          val,
          double          scaling,
          const long int* row_indices,
          const long int* col_indices,
          long int*       rows,
          long int*       cols,
          double*         vals,
          long int&       nnz)
      {
        const auto row_mapped = row_indices[static_cast<size_t>(row)];
        const auto col_mapped = col_indices[static_cast<size_t>(col)];
        if (col_mapped != INVALID_INDEX<long int>)
        {
          rows[static_cast<size_t>(nnz)] = row_mapped;
          cols[static_cast<size_t>(nnz)] = col_mapped;
          vals[static_cast<size_t>(nnz)] = scaling * val;
          nnz++;
        }
      }

      /**
       * @brief Enzyme sparse store
       *
       * @details This takes in a row, column and value and stores them in buffers
       *
       * @tparam ScalarT - scalar data type
       * @tparam IdxT - matrix index data type
       *
       * @param[in] val - value to be stored
       * @param[in] row - row to be stored
       * @param[in] col - column to be stored
       * @param[in] scaling - scaling factor for values
       * @param[in] res_indices - Global residual indices
       * @param[in] var_indices - Global variable indices
       * @param[in,out] rows - buffer where row will be stored
       * @param[in,out] cols - buffer where col will be stored
       * @param[in,out] vals - buffer where val will be stored
       * @param[in,out] nnz - number of nonzeros
       */
      template <typename ScalarT, typename IdxT>
      __attribute__((always_inline)) static void sparse_store(
          ScalarT     val,
          IdxT        row,
          IdxT        col,
          ScalarT     scaling,
          const IdxT* row_indices,
          const IdxT* col_indices,
          IdxT*       rows,
          IdxT*       cols,
          ScalarT*    vals,
          IdxT&       nnz)
      {
        if (val == 0.0)
          return;

        row /= sizeof(ScalarT);

        // this template nightmare is because __attribute__((enzyme_sparse_accumulate)) does not support templates yet
        if constexpr (std::is_same<IdxT, size_t>::value)
        {
          if constexpr (std::is_same<IdxT, float>::value)
            inner_store_float_size_t(row, col, val, scaling, row_indices, col_indices, rows, cols, vals, nnz);
          else
            inner_store_double_size_t(row, col, val, scaling, row_indices, col_indices, rows, cols, vals, nnz);
        }
        else if constexpr (std::is_same<IdxT, long int>::value)
        {
          if constexpr (std::is_same<IdxT, double>::value)
            inner_store_float_long_int(row, col, val, scaling, row_indices, col_indices, rows, cols, vals, nnz);
          else
            inner_store_double_long_int(row, col, val, scaling, row_indices, col_indices, rows, cols, vals, nnz);
        }
        else
        {
          assert(0 && "unsupported type");
        }
      }

      /**
       * @brief Enzyme sparse load
       *
       * @tparam ScalarT - scalar data type
       * @tparam IdxT - matrix index data type
       */
      template <typename ScalarT, typename IdxT>
      __attribute__((always_inline)) static ScalarT sparse_load(IdxT, IdxT, IdxT*, IdxT*, ScalarT*)
      {
        return 0.0;
      }

      /**
       * @brief Enzyme identity store
       *
       * @tparam ScalarT - scalar data type
       * @tparam IdxT - matrix index data type
       */
      template <typename ScalarT, typename IdxT>
      __attribute__((always_inline)) static void ident_store(ScalarT, IdxT, IdxT)
      {
        assert(0 && "should never store");
      }

      /**
       * @brief Enzyme identity load
       *
       * @tparam ScalarT - scalar data type
       * @tparam IdxT - matrix index data type
       */
      template <typename ScalarT, typename IdxT>
      __attribute__((always_inline)) static ScalarT ident_load(IdxT row, IdxT col)
      {
        row /= sizeof(ScalarT);
        return (ScalarT) (row == col);
      }

      /**
       * @brief Append an entry whose row and column have global indices
       *
       * @note __attribute__((enzyme_sparse_accumulate)) does not support templates yet
       *
       * @param[in] row - local row
       * @param[in] col - local column
       * @param[in] val - value to be stored
       * @param[in] row_indices - Global row of each local row
       * @param[in] col_indices - Global column of each local column
       * @param[in,out] entries - Entries in global indices
       */
      [[maybe_unused]] __attribute__((enzyme_sparse_accumulate)) static void inner_mapped_store_double_size_t(
          size_t        row,
          size_t        col,
          double        val,
          const size_t* row_indices,
          const size_t* col_indices,
          CooEntries*   entries)
      {
        const size_t row_mapped = row_indices[row];
        const size_t col_mapped = col_indices[col];
        if (row_mapped != INVALID_INDEX<size_t> && col_mapped != INVALID_INDEX<size_t>)
        {
          entries->rows.push_back(row_mapped);
          entries->cols.push_back(col_mapped);
          entries->values.push_back(val);
        }
      }

      /**
       * @brief Append a lower-triangle entry whose row and column have global indices
       *
       * @note __attribute__((enzyme_sparse_accumulate)) does not support templates yet
       *
       * @param[in] row - local row
       * @param[in] col - local column
       * @param[in] val - value to be stored
       * @param[in] row_indices - Global row of each local row
       * @param[in] col_indices - Global column of each local column
       * @param[in,out] entries - Entries in global indices
       */
      [[maybe_unused]] __attribute__((enzyme_sparse_accumulate)) static void inner_lower_store_double_size_t(
          size_t        row,
          size_t        col,
          double        val,
          const size_t* row_indices,
          const size_t* col_indices,
          CooEntries*   entries)
      {
        const size_t row_mapped = row_indices[row];
        const size_t col_mapped = col_indices[col];
        if (row_mapped != INVALID_INDEX<size_t> && col_mapped != INVALID_INDEX<size_t> && row_mapped >= col_mapped)
        {
          entries->rows.push_back(row_mapped);
          entries->cols.push_back(col_mapped);
          entries->values.push_back(val);
        }
      }

      /**
       * @brief Enzyme sparse store into `CooEntries`
       *
       * Drops rows and columns without a global index. The zero test marks
       * the store for Enzyme's sparsity analysis, which replaces the test by
       * the structural index set. Numerical zeros are therefore stored, and
       * the entries do not depend on values.
       *
       * @param[in] val - value to be stored
       * @param[in] row - row offset in bytes
       * @param[in] col - local column
       * @param[in] row_indices - Global row of each local row
       * @param[in] col_indices - Global column of each local column
       * @param[in,out] entries - Entries in global indices
       */
      [[maybe_unused]] __attribute__((always_inline)) static void mapped_store(double        val,
                                                                               size_t        row,
                                                                               size_t        col,
                                                                               const size_t* row_indices,
                                                                               const size_t* col_indices,
                                                                               CooEntries*   entries)
      {
        if (val == 0.0)
        {
          return;
        }

        row /= sizeof(double);
        inner_mapped_store_double_size_t(row, col, val, row_indices, col_indices, entries);
      }

      /**
       * @brief `mapped_store` restricted to the global lower triangle
       *
       * @param[in] val - value to be stored
       * @param[in] row - row offset in bytes
       * @param[in] col - local column
       * @param[in] row_indices - Global row of each local row
       * @param[in] col_indices - Global column of each local column
       * @param[in,out] entries - Entries in global indices
       */
      [[maybe_unused]] __attribute__((always_inline)) static void lower_store(double        val,
                                                                              size_t        row,
                                                                              size_t        col,
                                                                              const size_t* row_indices,
                                                                              const size_t* col_indices,
                                                                              CooEntries*   entries)
      {
        if (val == 0.0)
        {
          return;
        }

        row /= sizeof(double);
        inner_lower_store_double_size_t(row, col, val, row_indices, col_indices, entries);
      }

      /**
       * @brief Enzyme load for `mapped_store` and `lower_store` views
       *
       * Returns zero, so an accumulation stores its increment as a separate
       * entry.
       */
      [[maybe_unused]] __attribute__((always_inline)) static double mapped_load(size_t,
                                                                                size_t,
                                                                                const size_t*,
                                                                                const size_t*,
                                                                                CooEntries*)
      {
        return 0.0;
      }
    } // namespace Sparse
  } // namespace Enzyme
} // namespace GridKit
