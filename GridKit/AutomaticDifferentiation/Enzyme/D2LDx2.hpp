/**
 * @file D2LDx2.hpp
 * @brief Enzyme sparse Hessian of an optimization Lagrangian.
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
       * @brief Enzyme automatic differentiation Hessian evaluator: Lagrangian Hessian, d2L/dx2
       *
       * Hessian of \f$L = \sigma f(x) + \lambda^T g(x)\f$ by one forward sweep
       * per local variable over the reverse-mode gradient of \f$L\f$. Only
       * structural entries are stored, in global indices.
       *
       * @tparam ModelT - model type with internal and external local variables
       * and constraints
       */
      template <typename ModelT>
      struct D2LDx2
      {
        /**
         * @param[in] model - Pointer to the model to be differentiated
         * @param[in] n_var - Number of local variables
         * @param[in] variable_indices - Global index of each local variable
         * @param[in] x - Local variables
         * @param[in] sigma - Objective factor
         * @param[in] lambda - Multipliers of the local constraint contributions
         * @param[out] rows - Row indices
         * @param[out] cols - Column indices
         * @param[out] vals - Values
         * @param[out] nnz - Number of nonzeros
         */
        template <typename ScalarT, typename IdxT, typename RealT>
        static void eval(const ModelT*  model,
                         const size_t   n_var,
                         const IdxT*    variable_indices,
                         const ScalarT* x,
                         const RealT    sigma,
                         const RealT*   lambda,
                         IdxT*          rows,
                         IdxT*          cols,
                         ScalarT*       vals,
                         IdxT&          nnz)
        {
          std::array<ScalarT, ObjectiveT::NUM_INPUTS> gradient{};
          for (size_t var_i = 0; var_i < n_var; ++var_i)
          {
            // Sparse storage. @see LowerSparseStorage.hpp
            ScalarT* seed      = __enzyme_todense<ScalarT*>((void*) ident_load<ScalarT, IdxT>,
                                                       (void*) ident_store<ScalarT, IdxT>,
                                                       var_i);
            ScalarT* dgradient = __enzyme_todense<ScalarT*>((void*) sparse_load<ScalarT, IdxT>,
                                                            (void*) sparse_store<ScalarT, IdxT>,
                                                            var_i,
                                                            1.0, // value scaling
                                                            variable_indices,
                                                            variable_indices,
                                                            rows,
                                                            cols,
                                                            vals,
                                                            &nnz);

            __enzyme_fwddiff<void>((void*) lagrangianGradient<ScalarT, RealT>,
                                   enzyme_const,
                                   model,
                                   enzyme_dup,
                                   x,
                                   seed,
                                   enzyme_dupnoneed,
                                   gradient.data(),
                                   dgradient,
                                   enzyme_const,
                                   sigma,
                                   enzyme_const,
                                   lambda);
          }
        }

      private:
        using ObjectiveT   = ModelWrapper<ModelT, MemberFunctions::Objective>;
        using ConstraintsT = ModelWrapper<ModelT, MemberFunctions::Constraints>;

        template <typename ScalarT, typename RealT>
        [[gnu::always_inline]] inline static ScalarT lagrangian(const ModelT*  model,
                                                                const ScalarT* x,
                                                                const RealT    sigma,
                                                                const RealT*   lambda)
        {
          std::array<ScalarT, ConstraintsT::NUM_OUTPUTS> g{};
          ConstraintsT::eval(model, x, g.data());

          ScalarT value = sigma * ObjectiveT::eval(model, x);
          for (size_t i = 0; i < ConstraintsT::NUM_OUTPUTS; ++i)
          {
            value += lambda[i] * g[i];
          }
          return value;
        }

        /**
         * @brief Add the gradient of the Lagrangian to `gradient`
         *
         * Only the tangent of `gradient` is used, so it is not reset.
         */
        template <typename ScalarT, typename RealT>
        [[gnu::always_inline]] inline static void lagrangianGradient(const ModelT*  model,
                                                                     const ScalarT* x,
                                                                     ScalarT*       gradient,
                                                                     const RealT    sigma,
                                                                     const RealT*   lambda)
        {
          __enzyme_autodiff<void>((void*) lagrangian<ScalarT, RealT>,
                                  enzyme_const,
                                  model,
                                  enzyme_dup,
                                  x,
                                  gradient,
                                  enzyme_const,
                                  sigma,
                                  enzyme_const,
                                  lambda);
        }
      };
    } // namespace Sparse
  } // namespace Enzyme
} // namespace GridKit
