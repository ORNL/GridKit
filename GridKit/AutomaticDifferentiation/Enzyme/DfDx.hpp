/**
 * @file DfDx.hpp
 * @brief Enzyme gradient of an optimization objective.
 */

#pragma once

#include <algorithm>
#include <cstddef>

#include <GridKit/AutomaticDifferentiation/Enzyme/EnzymeDefinitions.hpp>
#include <GridKit/AutomaticDifferentiation/Enzyme/ModelWrappers.hpp>

namespace GridKit
{
  namespace Enzyme
  {
    namespace Sparse
    {
      /**
       * @brief Enzyme automatic differentiation gradient evaluator: objective gradient, df/dx
       *
       * One reverse sweep.
       *
       * @tparam ModelT - model type with `VARIABLE_SIZE` local variables
       */
      template <typename ModelT>
      struct DfDx
      {
        using ScalarT = typename ModelT::ScalarT;
        using RealT   = typename ModelT::RealT;

        /**
         * @param[in] model - Pointer to the model to be differentiated
         * @param[in] x - Local variables
         * @param[out] gradient - Gradient with respect to the local variables
         */
        static void eval(const ModelT* model, const ScalarT* x, RealT* gradient)
        {
          std::fill_n(gradient, ModelT::VARIABLE_SIZE, 0.0);

          __enzyme_autodiff<void>((void*) ModelWrapper<ModelT, MemberFunctions::Objective>::eval,
                                  enzyme_const,
                                  model,
                                  enzyme_dup,
                                  x,
                                  gradient);
        }
      };
    } // namespace Sparse
  } // namespace Enzyme
} // namespace GridKit
