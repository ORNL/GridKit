/**
 * @file D2LDx2.hpp
 * @brief Enzyme sparse Hessian of an optimization Lagrangian.
 */

#pragma once

#include <array>
#include <cstddef>
#include <type_traits>

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
       * Lower triangle of the Hessian of
       * \f$L = \sigma f(x) + \lambda^T g(x)\f$ by one forward sweep per local
       * variable over the reverse-mode gradient of \f$L\f$. Only structural
       * entries are stored, in global indices.
       *
       * @tparam ModelT - model type with `VARIABLE_SIZE` local variables and
       * `CONSTRAINT_SIZE` local constraints
       */
      template <typename ModelT>
      struct D2LDx2
      {
        using ScalarT = typename ModelT::ScalarT;
        using IdxT    = typename ModelT::IdxT;
        using RealT   = typename ModelT::RealT;

        static_assert(std::is_same_v<ScalarT, double> && std::is_same_v<IdxT, size_t>,
                      "D2LDx2 supports double values and size_t indices");

        /**
         * @param[in] model - Pointer to the model to be differentiated
         * @param[in] n_var - Number of local variables
         * @param[in] variable_indices - Global index of each local variable
         * @param[in] x - Local variables
         * @param[in] sigma - Objective factor
         * @param[in] lambda - Local constraint multipliers
         * @param[in,out] entries - Hessian entries
         */
        static void eval(const ModelT*  model,
                         const size_t   n_var,
                         const IdxT*    variable_indices,
                         const ScalarT* x,
                         const RealT    sigma,
                         const RealT*   lambda,
                         CooEntries*    entries)
        {
          std::array<ScalarT, ModelT::VARIABLE_SIZE> gradient{};
          for (size_t var_i = 0; var_i < n_var; ++var_i)
          {
            // Sparse storage. @see LowerSparseStorage.hpp
            ScalarT* seed      = __enzyme_todense<ScalarT*>((void*) ident_load<ScalarT, IdxT>,
                                                       (void*) ident_store<ScalarT, IdxT>,
                                                       var_i);
            ScalarT* dgradient = __enzyme_todense<ScalarT*>((void*) mapped_load,
                                                            (void*) lower_store,
                                                            var_i,
                                                            variable_indices,
                                                            variable_indices,
                                                            entries);

            __enzyme_fwddiff<void>((void*) lagrangianGradient,
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
        __attribute__((always_inline)) inline static ScalarT lagrangian(const ModelT*  model,
                                                                        const ScalarT* x,
                                                                        const RealT    sigma,
                                                                        const RealT*   lambda)
        {
          std::array<ScalarT, ModelT::CONSTRAINT_SIZE> g{};
          ModelWrapper<ModelT, MemberFunctions::Constraints>::eval(model, x, g.data());

          ScalarT value = sigma * ModelWrapper<ModelT, MemberFunctions::Objective>::eval(model, x);
          for (size_t i = 0; i < ModelT::CONSTRAINT_SIZE; ++i)
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
        __attribute__((always_inline)) inline static void lagrangianGradient(const ModelT*  model,
                                                                             const ScalarT* x,
                                                                             ScalarT*       gradient,
                                                                             const RealT    sigma,
                                                                             const RealT*   lambda)
        {
          __enzyme_autodiff<void>((void*) lagrangian,
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
