#pragma once

#include <vector>

/**
 * @brief Enzyme constants for activity analysis
 *
 */
extern int enzyme_dup;
extern int enzyme_const;
extern int enzyme_dupnoneed;

/**
 * @brief Residual wrapper around residual methods inside model classes
 *
 * @tparam ModelT - model type
 * @tparam ScalarT - scalar data type
 */
template <typename ModelT, typename ScalarT>
void residual_wrapper(ModelT* obj, ScalarT* y, ScalarT* f)
{
  obj->evaluateResidualLocally(y, f);
}

namespace GridKit
{
  namespace Enzyme
  {
    /**
     * @brief Enzyme fwddiff template
     *
     * @tparam T - return type
     * @tparam ModelT - model type
     */
    template <typename T, typename... ModelT>
    extern T __enzyme_fwddiff(void*, ModelT...) noexcept;

    /**
     * @brief Enzyme todense template
     *
     * @tparam T - return type
     * @tparam ModelT - model type
     */
    template <typename T, typename... ModelT>
    extern T __enzyme_todense(ModelT...) noexcept;

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
    __attribute__((enzyme_sparse_accumulate)) static void inner_storeflt(int64_t row, int64_t col, float val, std::vector<Triple<float>>& triplets)
    {
      triplets.emplace_back(row, col, val);
    }

    /**
     * @brief Enzyme sparse accumulation for double
     *
     */
    __attribute__((enzyme_sparse_accumulate)) static void inner_storedbl(int64_t row, int64_t col, double val, std::vector<Triple<double>>& triplets)
    {
      triplets.emplace_back(row, col, val);
    }

    /**
     * @brief Enzyme sparse store
     *
     * @tparam ScalarT - scalar data type
     */
    template <typename ScalarT>
    __attribute__((always_inline)) static void sparse_store(ScalarT val, int64_t idx, size_t i, std::vector<Triple<ScalarT>>& triplets)
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
    __attribute__((always_inline)) static ScalarT sparse_load(int64_t idx, size_t i, std::vector<Triple<ScalarT>>& triplets)
    {
      return 0.0;
    }

    /**
     * @brief Enzyme identity store
     *
     * @tparam ScalarT - scalar data type
     */
    template <typename ScalarT>
    __attribute__((always_inline)) static void ident_store(ScalarT, int64_t idx, size_t i)
    {
      assert(0 && "should never load");
    }

    /**
     * @brief Enzyme identity load
     *
     * @tparam ScalarT - scalar data type
     */
    template <typename ScalarT>
    __attribute__((always_inline)) static ScalarT ident_load(int64_t idx, size_t i)
    {
      idx /= sizeof(ScalarT);
      return (ScalarT) (idx == i);
    }

    /**
     * @brief Function that computes the Jacobian via automatic differentiation
     *
     * @tparam ModelT - model type
     * @tparam ScalarT - scalar data type
     * @tparam IdxT    - index data type
     */
    template <typename ModelT, class ScalarT, typename IdxT>
    __attribute__((noinline)) void EnzymeSparseModelJacobian(ModelT* model, IdxT n, ScalarT* input, GridKit::LinearAlgebra::COO_Matrix<ScalarT, IdxT>& jac)
    {
      std::vector<Triple<ScalarT>> triplets;
      for (IdxT i = 0; i < n; i++)
      {
        ScalarT* output   = __enzyme_todense<ScalarT*>((void*) ident_load<ScalarT>, (void*) ident_store<ScalarT>, i);
        ScalarT* d_output = __enzyme_todense<ScalarT*>((void*) sparse_load<ScalarT>, (void*) sparse_store<ScalarT>, i, &triplets);

        __enzyme_fwddiff<void>((void*) residual_wrapper<ModelT, ScalarT>,
                               enzyme_const,
                               model,
                               enzyme_dup,
                               input,
                               output,
                               enzyme_dupnoneed,
                               (ScalarT*) 0x1,
                               d_output);
      }

      std::vector<IdxT>    ctemp{};
      std::vector<IdxT>    rtemp{};
      std::vector<ScalarT> valtemp{};
      for (auto& tup : triplets)
      {
        rtemp.push_back(tup.row);
        ctemp.push_back(tup.col);
        valtemp.push_back(tup.val);
      }
      jac.setValues(rtemp, ctemp, valtemp);
    }
  } // namespace Enzyme
} // namespace GridKit
