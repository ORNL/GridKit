/**
 * @file SparseJacobian.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 *
 */

#pragma once

#include <GridKit/AutomaticDifferentiation/Enzyme/EnzymeDefinitions.hpp>
#include <GridKit/AutomaticDifferentiation/Enzyme/LowerSparseStorage.hpp>

namespace GridKit
{
  namespace Enzyme
  {
    namespace Sparse
    {
      /**
       * @brief Residual equation group keys
       *
       */
      enum class Equation
      {
        Internal,
        External
      };

      /**
       * @brief Differentiation variable group keys
       *
       */
      enum class Variable
      {
        Y,
        Yp,
        YExt,
        YpExt
      };

      /**
       * @brief Template definition for wrapper around residual methods inside model classes
       *
       * @tparam ModelT - model type
       * @tparam Equation - residual equation group key
       *
       */
      template <typename ModelT, Equation equation>
      struct ResidualWrapper
      {
      };

      /**
       * @brief Residual wrapper partial template specialization for the internal equations
       *
       */
      template <typename ModelT>
      struct ResidualWrapper<ModelT, Equation::Internal>
      {
        using ScalarT = typename ModelT::ScalarT;

        /**
         * @param[in] model - Pointer to the model to be differentiated
         * @param[in] y - Internal variables
         * @param[in] yp - Internal variable derivatives
         * @param[in] y_ext - External variables
         * @param[in] yp_ext - External variable derivatives
         * @param[out] f - Internal residual
         */
        static void eval(ModelT*        model,
                         const ScalarT* y,
                         const ScalarT* yp,
                         const ScalarT* y_ext,
                         const ScalarT* yp_ext,
                         ScalarT*       f)
        {
          model->evaluateInternalResidual(y, yp, y_ext, yp_ext, f);
        }
      };

      /**
       * @brief Residual wrapper partial template specialization for the external equations
       *
       */
      template <typename ModelT>
      struct ResidualWrapper<ModelT, Equation::External>
      {
        using ScalarT = typename ModelT::ScalarT;

        /**
         * @param[in] model - Pointer to the model to be differentiated
         * @param[in] y - Internal variables
         * @param[in] yp - Internal variable derivatives
         * @param[in] y_ext - External variables
         * @param[in] yp_ext - External variable derivatives
         * @param[out] f_ext - External residual
         */
        static void eval(ModelT*        model,
                         const ScalarT* y,
                         const ScalarT* yp,
                         const ScalarT* y_ext,
                         const ScalarT* yp_ext,
                         ScalarT*       f_ext)
        {
          model->evaluateExternalResidual(y, yp, y_ext, yp_ext, f_ext);
        }
      };

      /**
       * @brief Enzyme automatic differentiation Jacobian evaluator over one
       * equation group and one variable group
       *
       * @details The eight blocks of a model Jacobian are selected by the
       * (Equation, Variable) template key pair. Derivative variable groups
       * (Yp, YpExt) pass the integrator coefficient as the scaling so alpha
       * multiplies the derivative matrix exactly once.
       *
       * @tparam ModelT - model type
       * @tparam Equation - residual equation group key
       * @tparam Variable - differentiation variable group key
       */
      template <typename ModelT, Equation equation, Variable variable>
      struct SparseJacobian
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
         * @param[in] y_ext - External variables
         * @param[in] yp_ext - External variable derivatives
         * @param[out] rows - Row indices
         * @param[out] cols - Column indices
         * @param[out] vals - Values
         * @param[out] nnz - Number of nonzeros
         * @param[in] scaling - Value scaling, alpha for derivative variables
         */
        static void eval(ModelT*        model,
                         const size_t   n_res,
                         const size_t   n_var,
                         const IdxT*    res_indices,
                         const IdxT*    var_indices,
                         const ScalarT* y,
                         const ScalarT* yp,
                         const ScalarT* y_ext,
                         const ScalarT* yp_ext,
                         IdxT*          rows,
                         IdxT*          cols,
                         RealT*         vals,
                         IdxT&          nnz,
                         const RealT    scaling = 1.0)
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
                                                             scaling, // value scaling
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
              if constexpr (variable == Variable::Y)
              {
                __enzyme_fwddiff<void>((void*) ResidualWrapper<ModelT, equation>::eval,
                                       enzyme_const,
                                       model,
                                       enzyme_dup,
                                       y,
                                       output,
                                       enzyme_const,
                                       yp,
                                       enzyme_const,
                                       y_ext,
                                       enzyme_const,
                                       yp_ext,
                                       enzyme_dupnoneed,
                                       elementary_v.data(),
                                       d_output);
              }
              else if constexpr (variable == Variable::Yp)
              {
                __enzyme_fwddiff<void>((void*) ResidualWrapper<ModelT, equation>::eval,
                                       enzyme_const,
                                       model,
                                       enzyme_const,
                                       y,
                                       enzyme_dup,
                                       yp,
                                       output,
                                       enzyme_const,
                                       y_ext,
                                       enzyme_const,
                                       yp_ext,
                                       enzyme_dupnoneed,
                                       elementary_v.data(),
                                       d_output);
              }
              else if constexpr (variable == Variable::YExt)
              {
                __enzyme_fwddiff<void>((void*) ResidualWrapper<ModelT, equation>::eval,
                                       enzyme_const,
                                       model,
                                       enzyme_const,
                                       y,
                                       enzyme_const,
                                       yp,
                                       enzyme_dup,
                                       y_ext,
                                       output,
                                       enzyme_const,
                                       yp_ext,
                                       enzyme_dupnoneed,
                                       elementary_v.data(),
                                       d_output);
              }
              else
              {
                __enzyme_fwddiff<void>((void*) ResidualWrapper<ModelT, equation>::eval,
                                       enzyme_const,
                                       model,
                                       enzyme_const,
                                       y,
                                       enzyme_const,
                                       yp,
                                       enzyme_const,
                                       y_ext,
                                       enzyme_dup,
                                       yp_ext,
                                       output,
                                       enzyme_dupnoneed,
                                       elementary_v.data(),
                                       d_output);
              }
            }
          }
        }
      };
    } // namespace Sparse
  } // namespace Enzyme
} // namespace GridKit
