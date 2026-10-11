/**
 * @file Jacobian.hpp
 *
 */

#pragma once

#include <vector>

#include <GridKit/AutomaticDifferentiation/Enzyme/EnzymeDefinitions.hpp>
#include <GridKit/AutomaticDifferentiation/Enzyme/LowerSparseStorage.hpp>
#include <GridKit/Utilities/Enum.hpp>

namespace GridKit
{
  namespace Enzyme
  {
    namespace Sparse
    {
      /**
       * @brief Enzyme automatic differentiation of a component Jacobian.
       *
       * Computes dF/dy + alpha dF/dyp for F = [f | f_ext] over the columns
       * [y | y_ext]: f holds the rows the component owns, f_ext its
       * contributions to rows owned elsewhere. Own columns seed y and yp in
       * one pass, external columns seed y_ext.
       *
       * @note Enzyme rejects a call site with no stores, so a column loop is
       *       compiled only for a nonempty variable set.
       *
       * @tparam ModelT - model type with InternalVariablesT and ExternalVariablesT
       */
      template <typename ModelT>
      struct Jacobian
      {
        using ScalarT = typename ModelT::ScalarT;
        using IdxT    = typename ModelT::IdxT;
        using RealT   = typename ModelT::RealT;

        static constexpr size_t NUM_OWN = Utilities::enum_size<typename ModelT::InternalVariablesT>();
        static constexpr size_t NUM_EXT = Utilities::enum_size<typename ModelT::ExternalVariablesT>();

        /**
         * @brief Component residual F = [f | f_ext]
         *
         * @param[in] model - Pointer to the model to be differentiated
         * @param[in] y - Internal variables
         * @param[in] yp - Internal variable derivatives
         * @param[in] y_ext - External variables
         * @param[out] f - Internal residual
         * @param[out] f_ext - External residual
         */
        static void residual(ModelT*        model,
                             const ScalarT* y,
                             const ScalarT* yp,
                             const ScalarT* y_ext,
                             ScalarT*       f,
                             ScalarT*       f_ext)
        {
          model->evaluateInternalResidual(y, yp, y_ext, f);
          model->evaluateExternalResidual(y, yp, y_ext, f_ext);
        }

        /**
         * @param[in] model - Pointer to the model to be differentiated
         * @param[in] n - Number of internal variables and residuals
         * @param[in] n_ext - Number of external variables and residuals
         * @param[in] y - Internal variables
         * @param[in] yp - Internal variable derivatives
         * @param[in] y_ext - External variables
         * @param[in] alpha - Scaling of the derivative with respect to yp
         * @param[in] res_indices - Global residual indices of f
         * @param[in] var_indices - Global variable indices of y
         * @param[in] res_indices_ext - Global residual indices of f_ext
         * @param[in] var_indices_ext - Global variable indices of y_ext
         * @param[out] rows - Row indices
         * @param[out] cols - Column indices
         * @param[out] vals - Values
         * @param[out] nnz - Number of nonzeros
         */
        static void eval(ModelT*        model,
                         const size_t   n,
                         const size_t   n_ext,
                         const ScalarT* y,
                         const ScalarT* yp,
                         const ScalarT* y_ext,
                         const RealT    alpha,
                         const IdxT*    res_indices,
                         const IdxT*    var_indices,
                         const IdxT*    res_indices_ext,
                         const IdxT*    var_indices_ext,
                         IdxT*          rows,
                         IdxT*          cols,
                         RealT*         vals,
                         IdxT&          nnz)
        {
          // Only the tangents of the residuals are stored.
          std::vector<ScalarT> f(n);
          std::vector<ScalarT> f_ext(n_ext);

          if constexpr (NUM_OWN > 0)
          {
            for (size_t var_i = 0; var_i < n; ++var_i)
            {
              // Sparse storage. @see LowerSparseStorage.hpp
              ScalarT* d_y     = __enzyme_todense<ScalarT*>((void*) ident_load<ScalarT, IdxT>,
                                                        (void*) ident_store<ScalarT, IdxT>,
                                                        var_i);
              ScalarT* d_yp    = __enzyme_todense<ScalarT*>((void*) scaled_ident_load<ScalarT, IdxT>,
                                                         (void*) scaled_ident_store<ScalarT, IdxT>,
                                                         var_i,
                                                         alpha);
              ScalarT* d_f     = __enzyme_todense<ScalarT*>((void*) sparse_load<ScalarT, IdxT>,
                                                        (void*) sparse_store<ScalarT, IdxT>,
                                                        var_i,
                                                        res_indices,
                                                        var_indices,
                                                        rows,
                                                        cols,
                                                        vals,
                                                        &nnz);
              ScalarT* d_f_ext = __enzyme_todense<ScalarT*>((void*) sparse_load<ScalarT, IdxT>,
                                                            (void*) sparse_store<ScalarT, IdxT>,
                                                            var_i,
                                                            res_indices_ext,
                                                            var_indices,
                                                            rows,
                                                            cols,
                                                            vals,
                                                            &nnz);

              // Core automatic differentiation intrinsic that will be replaced by a derivative
              __enzyme_fwddiff<void>((void*) residual,
                                     enzyme_const,
                                     model,
                                     enzyme_dup,
                                     y,
                                     d_y,
                                     enzyme_dup,
                                     yp,
                                     d_yp,
                                     enzyme_const,
                                     y_ext,
                                     enzyme_dupnoneed,
                                     f.data(),
                                     d_f,
                                     enzyme_dupnoneed,
                                     f_ext.data(),
                                     d_f_ext);
            }
          }

          if constexpr (NUM_EXT > 0)
          {
            for (size_t var_i = 0; var_i < n_ext; ++var_i)
            {
              // Sparse storage. @see LowerSparseStorage.hpp
              ScalarT* d_y_ext = __enzyme_todense<ScalarT*>((void*) ident_load<ScalarT, IdxT>,
                                                            (void*) ident_store<ScalarT, IdxT>,
                                                            var_i);
              ScalarT* d_f     = __enzyme_todense<ScalarT*>((void*) sparse_load<ScalarT, IdxT>,
                                                        (void*) sparse_store<ScalarT, IdxT>,
                                                        var_i,
                                                        res_indices,
                                                        var_indices_ext,
                                                        rows,
                                                        cols,
                                                        vals,
                                                        &nnz);
              ScalarT* d_f_ext = __enzyme_todense<ScalarT*>((void*) sparse_load<ScalarT, IdxT>,
                                                            (void*) sparse_store<ScalarT, IdxT>,
                                                            var_i,
                                                            res_indices_ext,
                                                            var_indices_ext,
                                                            rows,
                                                            cols,
                                                            vals,
                                                            &nnz);

              // Core automatic differentiation intrinsic that will be replaced by a derivative
              __enzyme_fwddiff<void>((void*) residual,
                                     enzyme_const,
                                     model,
                                     enzyme_const,
                                     y,
                                     enzyme_const,
                                     yp,
                                     enzyme_dup,
                                     y_ext,
                                     d_y_ext,
                                     enzyme_dupnoneed,
                                     f.data(),
                                     d_f,
                                     enzyme_dupnoneed,
                                     f_ext.data(),
                                     d_f_ext);
            }
          }
        }
      };
    } // namespace Sparse
  } // namespace Enzyme
} // namespace GridKit
