/**
 * @file DgDx.hpp
 * @brief Enzyme sparse Jacobian of optimization constraints.
 */

#pragma once

#include <array>
#include <cstddef>

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
       * @brief Enzyme automatic differentiation Jacobian evaluator: constraint Jacobian, dg/dx
       *
       * One forward sweep per local variable. Only structural entries are
       * stored, in global indices.
       *
       * @tparam ModelT - model type with internal and external local variables
       * @tparam function - constraint functions, `Equalities` or `Inequalities`
       */
      template <typename ModelT, MemberFunctions function>
      struct DgDx
      {
        /**
         * @param[in] model - Pointer to the model to be differentiated
         * @param[in] n_var - Number of local variables
         * @param[in] row_indices - Global row of each local constraint function
         * @param[in] col_indices - Global column of each local variable
         * @param[in] x - Local variables
         * @param[out] rows - Row indices
         * @param[out] cols - Column indices
         * @param[out] vals - Values
         * @param[out] nnz - Number of nonzeros
         */
        template <typename ScalarT, typename IdxT>
        static void eval(const ModelT*  model,
                         const size_t   n_var,
                         const IdxT*    row_indices,
                         const IdxT*    col_indices,
                         const ScalarT* x,
                         IdxT*          rows,
                         IdxT*          cols,
                         ScalarT*       vals,
                         IdxT&          nnz)
        {
          using WrapperT = ModelWrapper<ModelT, function>;

          std::array<ScalarT, WrapperT::NUM_OUTPUTS> g{};
          for (size_t var_i = 0; var_i < n_var; ++var_i)
          {
            // Sparse storage. @see LowerSparseStorage.hpp
            ScalarT* seed = __enzyme_todense<ScalarT*>((void*) ident_load<ScalarT, IdxT>,
                                                       (void*) ident_store<ScalarT, IdxT>,
                                                       var_i);
            ScalarT* dg   = __enzyme_todense<ScalarT*>((void*) sparse_load<ScalarT, IdxT>,
                                                     (void*) sparse_store<ScalarT, IdxT>,
                                                     var_i,
                                                     1.0, // value scaling
                                                     row_indices,
                                                     col_indices,
                                                     rows,
                                                     cols,
                                                     vals,
                                                     &nnz);

            __enzyme_fwddiff<void>((void*) WrapperT::template eval<ScalarT>,
                                   enzyme_const,
                                   model,
                                   enzyme_dup,
                                   x,
                                   seed,
                                   enzyme_dupnoneed,
                                   g.data(),
                                   dg);
          }
        }
      };
    } // namespace Sparse
  } // namespace Enzyme
} // namespace GridKit
