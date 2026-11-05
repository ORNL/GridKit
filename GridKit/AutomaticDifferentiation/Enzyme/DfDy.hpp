/**
 * @file DfDy.hpp
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 *
 */

#pragma once

#include <array>
#include <cstddef>

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
       * @brief Enzyme automatic differentiation Jacobian evaluator: Internal Jacobian, df/dy
       *
       * @tparam ModelT - model type
       * @tparam MemberFunctions - member function parameter key
       */
      template <typename ModelT, MemberFunctions function>
      struct DfDy
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
                                     enzyme_dup,
                                     y,
                                     output, // Seed: column var_i of the identity
                                     enzyme_const,
                                     yp,
                                     enzyme_const,
                                     wb,
                                     enzyme_dupnoneed,
                                     (ScalarT*) 0x1, // Unneeded primal residual (dummy pointer)
                                     d_output);      // Column var_i of the Jacobian
            }
          }
        }

        /**
         * @param[in] model - Pointer to the model to be differentiated
         * @param[in] n_res - Number of residual functions
         * @param[in] n_var - Number of independent variables
         * @param[in] res_indices - Map from local residual indices to global indices
         * @param[in] var_indices - Map from local variable indices to global indices
         * @param[in] y - Internal variables
         * @param[in] yp - Internal variable derivatives
         * @param[in] wb - Bus variables
         * @param[in] ws - Signal variables
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
                         const ScalarT* ws,
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
                                     enzyme_dup,
                                     y,
                                     output, // Seed: column var_i of the identity
                                     enzyme_const,
                                     yp,
                                     enzyme_const,
                                     wb,
                                     enzyme_const,
                                     ws,
                                     enzyme_dupnoneed,
                                     (ScalarT*) 0x1, // Unneeded primal residual (dummy pointer)
                                     d_output);      // Column var_i of the Jacobian
            }
          }
        }
      };

      /**
       * @brief Enzyme automatic differentiation Jacobian evaluator: EMT residual, beta df/dy + alpha df/dyp
       *
       * Internal and external variables are seeded in separate sweeps.
       *
       * @tparam ModelT - Fixed-layout EMT model
       */
      template <typename ModelT>
      struct DfDy<ModelT, MemberFunctions::Residual>
      {
        /**
         * @param[in] model - Model to differentiate
         * @param[in] n_internal - Number of internal variables
         * @param[in] n_external - Number of external variables
         * @param[in] index - Global variable and residual of each local variable
         * @param[in] y - Local values
         * @param[in] yp - Local derivatives
         * @param[in] beta - Scale of df/dy
         * @param[in] alpha - Scale of df/dyp
         * @param[out] rows - Row indices
         * @param[out] cols - Column indices
         * @param[out] vals - Values
         * @param[out] nnz - Number of entries
         */
        template <typename ScalarT, typename IdxT>
        static void eval(const ModelT*  model,
                         const size_t   n_internal,
                         const size_t   n_external,
                         const IdxT*    index,
                         const ScalarT* y,
                         const ScalarT* yp,
                         const ScalarT  beta,
                         const ScalarT  alpha,
                         IdxT*          rows,
                         IdxT*          cols,
                         ScalarT*       vals,
                         IdxT&          nnz)
        {
          using WrapperT = ModelWrapper<ModelT, MemberFunctions::Residual>;

          std::array<ScalarT, WrapperT::NUM_INTERNAL + WrapperT::NUM_EXTERNAL> f{};

          const ScalarT* w  = y + WrapperT::NUM_INTERNAL;
          const ScalarT* wp = yp + WrapperT::NUM_INTERNAL;

          if constexpr (WrapperT::NUM_INTERNAL > 0)
          {
            for (size_t j = 0; j < n_internal; ++j)
            {
              ScalarT* dy  = __enzyme_todense<ScalarT*>((void*) scaled_ident_load<ScalarT, IdxT>,
                                                       (void*) scaled_ident_store<ScalarT, IdxT>,
                                                       j,
                                                       beta);
              ScalarT* dyp = __enzyme_todense<ScalarT*>((void*) scaled_ident_load<ScalarT, IdxT>,
                                                        (void*) scaled_ident_store<ScalarT, IdxT>,
                                                        j,
                                                        alpha);
              ScalarT* df  = __enzyme_todense<ScalarT*>((void*) sparse_load<ScalarT, IdxT>,
                                                       (void*) sparse_store<ScalarT, IdxT>,
                                                       j,
                                                       1.0,
                                                       index,
                                                       index,
                                                       rows,
                                                       cols,
                                                       vals,
                                                       &nnz);

              __enzyme_fwddiff<void>((void*) WrapperT::template eval<ScalarT>,
                                     enzyme_const,
                                     model,
                                     enzyme_dup,
                                     y,
                                     dy,
                                     enzyme_dup,
                                     yp,
                                     dyp,
                                     enzyme_const,
                                     w,
                                     enzyme_const,
                                     wp,
                                     enzyme_dupnoneed,
                                     f.data(),
                                     df);
            }
          }

          if constexpr (WrapperT::NUM_EXTERNAL > 0)
          {
            for (size_t j = 0; j < n_external; ++j)
            {
              ScalarT* dw  = __enzyme_todense<ScalarT*>((void*) scaled_ident_load<ScalarT, IdxT>,
                                                       (void*) scaled_ident_store<ScalarT, IdxT>,
                                                       j,
                                                       beta);
              ScalarT* dwp = __enzyme_todense<ScalarT*>((void*) scaled_ident_load<ScalarT, IdxT>,
                                                        (void*) scaled_ident_store<ScalarT, IdxT>,
                                                        j,
                                                        alpha);
              ScalarT* df  = __enzyme_todense<ScalarT*>((void*) sparse_load<ScalarT, IdxT>,
                                                       (void*) sparse_store<ScalarT, IdxT>,
                                                       j,
                                                       1.0,
                                                       index,
                                                       index + WrapperT::NUM_INTERNAL,
                                                       rows,
                                                       cols,
                                                       vals,
                                                       &nnz);

              __enzyme_fwddiff<void>((void*) WrapperT::template eval<ScalarT>,
                                     enzyme_const,
                                     model,
                                     enzyme_const,
                                     y,
                                     enzyme_const,
                                     yp,
                                     enzyme_dup,
                                     w,
                                     dw,
                                     enzyme_dup,
                                     wp,
                                     dwp,
                                     enzyme_dupnoneed,
                                     f.data(),
                                     df);
            }
          }
        }
      };
    } // namespace Sparse
  } // namespace Enzyme
} // namespace GridKit
