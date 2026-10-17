/**
 * @file DhDwb.hpp
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 *
 */

#pragma once

#include <GridKit/AutomaticDifferentiation/Enzyme/EnzymeDefinitions.hpp>
#include <GridKit/AutomaticDifferentiation/Enzyme/LowerSparseStorage.hpp>
#include <GridKit/AutomaticDifferentiation/Enzyme/ModelWrappers.hpp>
#include <GridKit/Model/EvaluatorTraits.hpp>

namespace GridKit
{
  namespace Enzyme
  {
    namespace Sparse
    {
      /**
       * @brief Enzyme automatic differentiation Jacobian evaluator: dh/dwb
       *
       * @tparam ModelT - model type
       * @tparam MemberFunctions - member function parameter key
       */
      template <typename ModelT, MemberFunctions function>
      struct DhDwb
      {
        using ScalarT = typename GridKit::Model::EvaluatorTraits<ModelT>::ScalarT;
        using IdxT    = typename GridKit::Model::EvaluatorTraits<ModelT>::IdxT;
        using RealT   = typename GridKit::Model::EvaluatorTraits<ModelT>::RealT;

        /**
         * @param[in] model - Pointer to the model to be differentiated
         * @param[in] n_res - Number of residual functions
         * @param[in] n_var - Number of independent variables
         * @param[in] res_indices - Global residual indices
         * @param[in] var_indices - Global variable indices
         * @param[in] y - Internal variables
         * @param[in] yp - Internal variable derivatives
         * @param[in] wb - Bus variables
         * @param[out] rows - Row indices
         * @param[out] cols - Column indices
         * @param[out] vals - Values
         * @param[out] nnz - Number of nonzeros
         */
        static void eval(ModelT*        model,
                         const size_t   n_res,
                         const size_t   n_var,
                         const IdxT*    res_indices,
                         const IdxT*    var_indices,
                         const ScalarT* y,
                         const ScalarT* yp,
                         const ScalarT* wb,
                         IdxT*          rows,
                         IdxT*          cols,
                         RealT*         vals,
                         IdxT&          nnz)
        {
          if (n_res > 0 && n_var > 0)
          {
            for (size_t var_i = 0; var_i < n_var; ++var_i)
            {
              // Sparse storage. @see LowerSparseStorage.hpp
              ScalarT* output   = __enzyme_todense<ScalarT*>((void*) ident_load<ScalarT, IdxT>,
                                                           (void*) ident_store<ScalarT, IdxT>,
                                                           var_i);
              ScalarT* d_output = __enzyme_todense<ScalarT*>((void*) sparse_load<ScalarT, IdxT>,
                                                             (void*) sparse_store<ScalarT, IdxT>,
                                                             var_i,
                                                             1.0, // value scaling
                                                             res_indices,
                                                             var_indices,
                                                             rows,
                                                             cols,
                                                             vals,
                                                             &nnz);

              // Core automatic differentiation intrinsic that will be replaced by a derivative
              // The dummy primal pointer is safe only while residuals never read back their output
              __enzyme_fwddiff<void>((void*) ModelWrapper<ModelT, function>::eval,
                                     enzyme_const,
                                     model,
                                     enzyme_const,
                                     y,
                                     enzyme_const,
                                     yp,
                                     enzyme_dup,
                                     wb,
                                     output, // Seed: column var_i of the identity
                                     enzyme_dupnoneed,
                                     (ScalarT*) 0x1, // Unneeded primal residual (dummy pointer)
                                     d_output);      // Column var_i of the Jacobian
            }
          }
        }
      };
    } // namespace Sparse
  } // namespace Enzyme
} // namespace GridKit
