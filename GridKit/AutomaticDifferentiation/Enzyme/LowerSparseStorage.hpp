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
       * @tparam scalar_type - scalar data type
       * @tparam index_type - matrix index data type
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
      template <typename scalar_type, typename index_type>
      __attribute__((always_inline)) static void sparse_store(
          scalar_type       val,
          index_type        row,
          index_type        col,
          scalar_type       scaling,
          const index_type* row_indices,
          const index_type* col_indices,
          index_type*       rows,
          index_type*       cols,
          scalar_type*      vals,
          index_type&       nnz)
      {
        if (val == 0.0)
          return;

        row /= sizeof(scalar_type);

        // this template nightmare is because __attribute__((enzyme_sparse_accumulate)) does not support templates yet
        if constexpr (std::is_same<index_type, size_t>::value)
        {
          if constexpr (std::is_same<index_type, float>::value)
            inner_store_float_size_t(row, col, val, scaling, row_indices, col_indices, rows, cols, vals, nnz);
          else
            inner_store_double_size_t(row, col, val, scaling, row_indices, col_indices, rows, cols, vals, nnz);
        }
        else if constexpr (std::is_same<index_type, long int>::value)
        {
          if constexpr (std::is_same<index_type, double>::value)
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
       * @tparam scalar_type - scalar data type
       * @tparam index_type - matrix index data type
       */
      template <typename scalar_type, typename index_type>
      __attribute__((always_inline)) static scalar_type sparse_load(index_type, index_type, index_type*, index_type*, scalar_type*)
      {
        return 0.0;
      }

      /**
       * @brief Enzyme identity store
       *
       * @tparam scalar_type - scalar data type
       * @tparam index_type - matrix index data type
       */
      template <typename scalar_type, typename index_type>
      __attribute__((always_inline)) static void ident_store(scalar_type, index_type, index_type)
      {
        assert(0 && "should never store");
      }

      /**
       * @brief Enzyme identity load
       *
       * @tparam scalar_type - scalar data type
       * @tparam index_type - matrix index data type
       */
      template <typename scalar_type, typename index_type>
      __attribute__((always_inline)) static scalar_type ident_load(index_type row, index_type col)
      {
        row /= sizeof(scalar_type);
        return (scalar_type) (row == col);
      }
    } // namespace Sparse
  } // namespace Enzyme
} // namespace GridKit
