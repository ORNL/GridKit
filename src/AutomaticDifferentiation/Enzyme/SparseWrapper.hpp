#pragma once

#include <vector>

/**
 * @brief Enzyme constants for activity analysis
 *
 */
extern int enzyme_dup;
extern int enzyme_const;
extern int enzyme_dupnoneed;

namespace GridKit
{
  namespace Enzyme
  {
    namespace Sparse
    {
      /**
       * @brief Residual wrapper around residual methods inside model classes
       *
       * @tparam ModelT - model type
       * @tparam ScalarT - scalar data type
       */
      template <typename ModelT, typename ScalarT>
      void residual_wrapper(ModelT* obj, ScalarT* y, ScalarT* yp, ScalarT* f)
      {
        obj->evaluateResidualLocally(y, yp, f);
      }

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

      /**
       * @brief Function that computes the Jacobian via automatic differentiation
       *
       * @tparam ModelT - model type
       * @tparam ScalarT - scalar data type
       * @tparam IdxT    - matrix index data type
       */
      template <typename ModelT, class ScalarT, typename IdxT>
      void ModelJacobian(ModelT* model, size_t n_res, size_t n_state, ScalarT* y, ScalarT* yp, GridKit::LinearAlgebra::COO_Matrix<ScalarT, IdxT>& jac)
      {
        std::vector<Triple<ScalarT>> triplets;
        std::vector<ScalarT>         elementary_v(n_state);
        for (size_t res_i = 0; res_i < n_res; ++res_i)
        {
          // Sparse storage
          ScalarT* output   = __enzyme_todense<ScalarT*>((void*) ident_load<ScalarT>, (void*) ident_store<ScalarT>, res_i);
          ScalarT* d_output = __enzyme_todense<ScalarT*>((void*) sparse_load<ScalarT>, (void*) sparse_store<ScalarT>, res_i, &triplets);

          // Elementary vector for Jacobian-vector product
          for (size_t state_i = 0; state_i < n_state; ++state_i)
          {
            elementary_v[state_i] = 0.0;
          }
          elementary_v[res_i] = 1.0;

          // Autodiff
          __enzyme_fwddiff<void>((void*) residual_wrapper<ModelT, ScalarT>,
                                 enzyme_const,
                                 model,
                                 enzyme_dup,
                                 y,
                                 output,
                                 enzyme_const,
                                 yp,
                                 enzyme_dupnoneed,
                                 elementary_v.data(),
                                 d_output);
        }

        // Store result
        std::vector<IdxT>    ctemp{};
        std::vector<IdxT>    rtemp{};
        std::vector<ScalarT> valtemp{};
        for (auto& tup : triplets)
        {
          rtemp.push_back(static_cast<IdxT>(tup.row));
          ctemp.push_back(static_cast<IdxT>(tup.col));
          valtemp.push_back(tup.val);
        }
        jac.setValues(rtemp, ctemp, valtemp);
      }
    } // namespace Sparse
  } // namespace Enzyme
} // namespace GridKit
