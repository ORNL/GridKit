/**
 * @file DhDy.hpp
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 *
 */

#pragma once

#include <GridKit/AutomaticDifferentiation/Enzyme/EnzymeDefinitions.hpp>
#include <GridKit/AutomaticDifferentiation/Enzyme/LowerSparseStorage.hpp>
#include <GridKit/AutomaticDifferentiation/Enzyme/ModelWrappers.hpp>

namespace GridKit
{
  namespace Enzyme
  {
    namespace Sparse
    {
      /**
       * @brief Enzyme automatic differentiation Jacobian evaluator: dh/dy
       *
       * @tparam ModelT - model type
       * @tparam MemberFunctions - member function parameter key
       */
      template <typename ModelT, MemberFunctions function>
      struct DhDy
      {
        using ScalarT = typename ModelT::ScalarT;
        using IdxT    = typename ModelT::IdxT;
        using RealT   = typename ModelT::RealT;

        /**
         * @param[in] model - Pointer to the model to be differentiated
         * @param[in] n_res - Number of residual functions
         * @param[in] n_var - Number of independent variables
         * @param[in] res_indices - Global residual indices
         * @param[in] var_indices - Global variable indices
         * @param[in] y - Internal variables
         * @param[in] yp - Internal variable derivatives
         * @param[in] wb - Bus variables
         * @param[out] nnz - Number of nonzeros
         */
        static void eval(ModelT*     model,
                         size_t      n_res,
                         size_t      n_var,
                         const IdxT* res_indices,
                         const IdxT* var_indices,
                         ScalarT*    y,
                         ScalarT*    yp,
                         ScalarT*    wb,
                         IdxT*       rows,
                         IdxT*       cols,
                         RealT*      vals,
                         IdxT&       nnz)
        {
          if (n_res > 0 && n_var > 0)
          {
            std::vector<ScalarT> elementary_v(n_var);
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

              // Elementary vector for Jacobian-vector product
              std::ranges::fill(elementary_v, 0.0);
              elementary_v[var_i] = 1.0;

              // Core automatic differentiation intrinsic that will be replaced by a derivative
              __enzyme_fwddiff<void>((void*) ModelWrapper<ModelT, function>::eval,
                                     enzyme_const,
                                     model,
                                     enzyme_dup,
                                     y,
                                     output,
                                     enzyme_const,
                                     yp,
                                     enzyme_const,
                                     wb,
                                     enzyme_dupnoneed,
                                     elementary_v.data(),
                                     d_output);
            }
          }
        }
      };
    } // namespace Sparse
  } // namespace Enzyme
} // namespace GridKit
