/**
 * @file DfDy.hpp
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 *
 */

#pragma once

#include <GridKit/AutomaticDifferentiation/Enzyme/EnzymeDefinitions.hpp>
#include <GridKit/AutomaticDifferentiation/Enzyme/LowerSparseStorage.hpp>
#include <GridKit/AutomaticDifferentiation/Enzyme/ModelWrappers.hpp>
#include <GridKit/LinearAlgebra/SparseMatrix/COO_Matrix.hpp>
#include <GridKit/ScalarTraits.hpp>

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
       * @tparam ScalarT - scalar data type
       * @tparam IdxT - matrix index data type
       */
      template <typename ModelT, MemberFunctions function, class ScalarT, typename IdxT>
      struct DfDy
      {
        using RealT   = typename GridKit::ScalarTraits<ScalarT>::RealT;
        using MatrixT = GridKit::LinearAlgebra::COO_Matrix<RealT, IdxT>;

        /**
         * @param[in] model - Pointer to the model to be differentiated
         * @param[in] n_res - Number of residual functions
         * @param[in] n_var - Number of independent variables
         * @param[in] res_indices - Global residual indices
         * @param[in] var_indices - Global variable indices
         * @param[in] y - Internal variables
         * @param[in] yp - Internal variable derivatives
         * @param[in] wb - Bus variables
         * @param[in,out] jac - Jacobian
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
                         MatrixT&    jac)
        {
          if (n_res > 0 && n_var > 0)
          {
            // df/dy
            std::vector<ScalarT> elementary_v(n_var);
            IdxT                 nnz = 0;
            for (size_t var_i = 0; var_i < n_var; ++var_i)
            {
              // Sparse storage. @see LowerSparseStorage.hpp
              ScalarT* output   = __enzyme_todense<ScalarT*>((void*) ident_load<ScalarT, IdxT>,
                                                           (void*) ident_store<ScalarT, IdxT>,
                                                           var_i);
              ScalarT* d_output = __enzyme_todense<ScalarT*>((void*) sparse_load<ScalarT, IdxT>,
                                                             (void*) sparse_store<ScalarT, IdxT>,
                                                             var_i,
                                                             res_indices,
                                                             var_indices,
                                                             rows,
                                                             cols,
                                                             vals,
                                                             &nnz);

              // Elementary vector for Jacobian-vector product
              std::ranges::fill(elementary_v, 0.0);
              elementary_v[var_i] = 1.0;

              // Core automatic differentiaation intrinsic that will be replaced by a derivative
              __enzyme_fwddiff<void>((void*) ModelWrapper<ModelT, function, ScalarT>::eval,
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

            // Store result
            jac.setValues(1.0, rows, cols, vals, nnz); //< @todo: Update once sparse storage format changes

            // There is no df/dy' when alpha is not passed as an argument
            // @todo: Implement a generic way to identify these cases at compile time
          }
        }

        /**
         * @param[in] model - Pointer to the model to be differentiated
         * @param[in] n_res - Number of residual functions
         * @param[in] n_var - Number of independent variables
         * @param[in] res_indices - Global residual indices
         * @param[in] var_indices - Global variable indices
         * @param[in] y - Internal variables
         * @param[in] yp - Internal variable derivatives
         * @param[in] wb - Bus variables
         * @param[in] alpha - Time derivative jacobian coefficient
         * @param[in,out] jac - Jacobian
         */
        static void eval(ModelT*     model,
                         size_t      n_res,
                         size_t      n_var,
                         const IdxT* res_indices,
                         const IdxT* var_indices,
                         ScalarT*    y,
                         ScalarT*    yp,
                         ScalarT*    wb,
                         RealT       alpha,
                         IdxT*       rows,
                         IdxT*       cols,
                         RealT*      vals,
                         MatrixT&    jac)
        {
          if (n_res > 0 && n_var > 0)
          {
            // df/dy
            std::vector<ScalarT> elementary_v(n_var);
            IdxT                 nnz = 0;
            for (size_t var_i = 0; var_i < n_var; ++var_i)
            {
              // Sparse storage. @see LowerSparseStorage.hpp
              ScalarT* output   = __enzyme_todense<ScalarT*>((void*) ident_load<ScalarT, IdxT>,
                                                           (void*) ident_store<ScalarT, IdxT>,
                                                           var_i);
              ScalarT* d_output = __enzyme_todense<ScalarT*>((void*) sparse_load<ScalarT, IdxT>,
                                                             (void*) sparse_store<ScalarT, IdxT>,
                                                             var_i,
                                                             res_indices,
                                                             var_indices,
                                                             rows,
                                                             cols,
                                                             vals,
                                                             &nnz);

              // Elementary vector for Jacobian-vector product
              std::ranges::fill(elementary_v, 0.0);
              elementary_v[var_i] = 1.0;

              // Core automatic differentiaation intrinsic that will be replaced by a derivative
              __enzyme_fwddiff<void>((void*) ModelWrapper<ModelT, function, ScalarT>::eval,
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

            // Store result
            jac.setValues(1.0, rows, cols, vals, nnz); //< @todo: Update once sparse storage format changes

            // df/dy'
            nnz = 0;
            for (size_t var_i = 0; var_i < n_var; ++var_i)
            {
              // Sparse storage. @see LowerSparseStorage.hpp
              ScalarT* output   = __enzyme_todense<ScalarT*>((void*) ident_load<ScalarT, IdxT>,
                                                           (void*) ident_store<ScalarT, IdxT>,
                                                           var_i);
              ScalarT* d_output = __enzyme_todense<ScalarT*>((void*) sparse_load<ScalarT, IdxT>,
                                                             (void*) sparse_store<ScalarT, IdxT>,
                                                             var_i,
                                                             res_indices,
                                                             var_indices,
                                                             rows,
                                                             cols,
                                                             vals,
                                                             &nnz);

              // Elementary vector for Jacobian-vector product
              std::ranges::fill(elementary_v, 0.0);
              elementary_v[var_i] = 1.0;

              // Core automatic differentiaation intrinsic that will be replaced by a derivative
              __enzyme_fwddiff<void>((void*) ModelWrapper<ModelT, function, ScalarT>::eval,
                                     enzyme_const,
                                     model,
                                     enzyme_const,
                                     y,
                                     enzyme_dup,
                                     yp,
                                     output,
                                     enzyme_const,
                                     wb,
                                     enzyme_dupnoneed,
                                     elementary_v.data(),
                                     d_output);
            }

            // Store result
            jac.setValues(alpha, rows, cols, vals, nnz); //< @todo: Update once sparse storage format changes
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
         * @param[in] alpha - Time derivative jacobian coefficient
         * @param[in,out] jac - Jacobian
         */
        static void eval(ModelT*     model,
                         size_t      n_res,
                         size_t      n_var,
                         const IdxT* res_indices,
                         const IdxT* var_indices,
                         ScalarT*    y,
                         ScalarT*    yp,
                         ScalarT*    wb,
                         ScalarT*    ws,
                         RealT       alpha,
                         IdxT*       rows,
                         IdxT*       cols,
                         RealT*      vals,
                         MatrixT&    jac)
        {
          if (n_res > 0 && n_var > 0)
          {
            // df/dy
            std::vector<ScalarT> elementary_v(n_var);
            IdxT                 nnz = 0;
            for (size_t var_i = 0; var_i < n_var; ++var_i)
            {
              // Sparse storage. @see LowerSparseStorage.hpp
              ScalarT* output   = __enzyme_todense<ScalarT*>((void*) ident_load<ScalarT, IdxT>,
                                                           (void*) ident_store<ScalarT, IdxT>,
                                                           var_i);
              ScalarT* d_output = __enzyme_todense<ScalarT*>((void*) sparse_load<ScalarT, IdxT>,
                                                             (void*) sparse_store<ScalarT, IdxT>,
                                                             var_i,
                                                             res_indices,
                                                             var_indices,
                                                             rows,
                                                             cols,
                                                             vals,
                                                             &nnz);

              // Elementary vector for Jacobian-vector product
              std::ranges::fill(elementary_v, 0.0);
              elementary_v[var_i] = 1.0;

              // Core automatic differentiaation intrinsic that will be replaced by a derivative
              __enzyme_fwddiff<void>((void*) ModelWrapper<ModelT, function, ScalarT>::eval,
                                     enzyme_const,
                                     model,
                                     enzyme_dup,
                                     y,
                                     output,
                                     enzyme_const,
                                     yp,
                                     enzyme_const,
                                     wb,
                                     enzyme_const,
                                     ws,
                                     enzyme_dupnoneed,
                                     elementary_v.data(),
                                     d_output);
            }

            // Store result
            jac.setValues(1.0, rows, cols, vals, nnz); //< @todo: Update once sparse storage format changes

            // df/dy'
            nnz = 0;
            for (size_t var_i = 0; var_i < n_var; ++var_i)
            {
              // Sparse storage. @see LowerSparseStorage.hpp
              ScalarT* output   = __enzyme_todense<ScalarT*>((void*) ident_load<ScalarT, IdxT>,
                                                           (void*) ident_store<ScalarT, IdxT>,
                                                           var_i);
              ScalarT* d_output = __enzyme_todense<ScalarT*>((void*) sparse_load<ScalarT, IdxT>,
                                                             (void*) sparse_store<ScalarT, IdxT>,
                                                             var_i,
                                                             res_indices,
                                                             var_indices,
                                                             rows,
                                                             cols,
                                                             vals,
                                                             &nnz);

              // Elementary vector for Jacobian-vector product
              std::ranges::fill(elementary_v, 0.0);
              elementary_v[var_i] = 1.0;

              // Core automatic differentiaation intrinsic that will be replaced by a derivative
              __enzyme_fwddiff<void>((void*) ModelWrapper<ModelT, function, ScalarT>::eval,
                                     enzyme_const,
                                     model,
                                     enzyme_const,
                                     y,
                                     enzyme_dup,
                                     yp,
                                     output,
                                     enzyme_const,
                                     wb,
                                     enzyme_const,
                                     ws,
                                     enzyme_dupnoneed,
                                     elementary_v.data(),
                                     d_output);
            }

            // Store result
            jac.setValues(alpha, rows, cols, vals, nnz); //< @todo: Update once sparse storage format changes
          }
        }
      };
    } // namespace Sparse
  } // namespace Enzyme
} // namespace GridKit
