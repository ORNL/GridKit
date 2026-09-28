/**
 * @file DgDx.hpp
 * @brief Enzyme sparse Jacobian of optimization constraints.
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
       * @brief Enzyme automatic differentiation Jacobian evaluator: constraint Jacobian, dg/dx
       *
       * One forward sweep per local variable. Only structural entries are
       * stored, in global indices.
       *
       * @tparam ModelT - model type with `CONSTRAINT_SIZE` local constraints
       */
      template <typename ModelT>
      struct DgDx
      {
        using ScalarT = typename ModelT::ScalarT;
        using IdxT    = typename ModelT::IdxT;

        static_assert(std::is_same_v<ScalarT, double> && std::is_same_v<IdxT, size_t>,
                      "DgDx supports double values and size_t indices");

        /**
         * @param[in] model - Pointer to the model to be differentiated
         * @param[in] n_var - Number of local variables
         * @param[in] constraint_indices - Global row of each local constraint
         * @param[in] variable_indices - Global column of each local variable
         * @param[in] x - Local variables
         * @param[in,out] entries - Jacobian entries
         */
        static void eval(const ModelT*  model,
                         const size_t   n_var,
                         const IdxT*    constraint_indices,
                         const IdxT*    variable_indices,
                         const ScalarT* x,
                         CooEntries*    entries)
        {
          std::array<ScalarT, ModelT::CONSTRAINT_SIZE> g{};
          for (size_t var_i = 0; var_i < n_var; ++var_i)
          {
            // Sparse storage. @see LowerSparseStorage.hpp
            ScalarT* seed = __enzyme_todense<ScalarT*>((void*) ident_load<ScalarT, IdxT>,
                                                       (void*) ident_store<ScalarT, IdxT>,
                                                       var_i);
            ScalarT* dg   = __enzyme_todense<ScalarT*>((void*) mapped_load,
                                                     (void*) mapped_store,
                                                     var_i,
                                                     constraint_indices,
                                                     variable_indices,
                                                     entries);

            __enzyme_fwddiff<void>((void*) ModelWrapper<ModelT, MemberFunctions::Constraints>::eval,
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
