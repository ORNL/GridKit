/**
 * @file LowerSparseStorage.hpp
 *
 * @details This file contains functions used by Enzyme to store sparse Jacobians.
 *
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 *
 */

#pragma once

#include <vector>

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
          vals[static_cast<size_t>(nnz)] = val;
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
          vals[static_cast<size_t>(nnz)] = val;
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
          vals[static_cast<size_t>(nnz)] = val;
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
          vals[static_cast<size_t>(nnz)] = val;
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
       * @param[in] res_indices - Global residual indices
       * @param[in] var_indices - Global variable indices
       * @param[in,out] rows - buffer where row will be stored
       * @param[in,out] cols - buffer where col will be stored
       * @param[in,out] vals - buffer where val will be stored
       * @param[in,out] nnz - number of nonzeros
       */
      template <typename ScalarT, typename IdxT>
      FORCE_INLINE static void sparse_store(
          ScalarT     val,
          IdxT        row,
          IdxT        col,
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
            inner_store_float_size_t(row, col, val, row_indices, col_indices, rows, cols, vals, nnz);
          else
            inner_store_double_size_t(row, col, val, row_indices, col_indices, rows, cols, vals, nnz);
        }
        else if constexpr (std::is_same<IdxT, long int>::value)
        {
          if constexpr (std::is_same<IdxT, double>::value)
            inner_store_float_long_int(row, col, val, row_indices, col_indices, rows, cols, vals, nnz);
          else
            inner_store_double_long_int(row, col, val, row_indices, col_indices, rows, cols, vals, nnz);
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
      FORCE_INLINE static ScalarT sparse_load(IdxT, IdxT, IdxT*, IdxT*, ScalarT*)
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
      FORCE_INLINE static void ident_store(ScalarT, IdxT, IdxT)
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
      FORCE_INLINE static ScalarT ident_load(IdxT row, IdxT col)
      {
        row /= sizeof(ScalarT);
        return (ScalarT) (row == col);
      }
    } // namespace Sparse
  } // namespace Enzyme
} // namespace GridKit
