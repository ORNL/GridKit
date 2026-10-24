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
       * @brief Enzyme objective gradient, df/dx, by one reverse sweep
       *
       * @tparam ModelT - optimization component model
       */
      template <typename ModelT>
      struct DfDx
      {
        /**
         * @param[in] model - Pointer to the model to be differentiated
         * @param[in] x - Local variables
         * @param[out] gradient - Gradient with respect to the local variables
         */
        template <typename ScalarT, typename RealT>
        static void eval(const ModelT* model, const ScalarT* x, RealT* gradient)
        {
          using WrapperT = ModelWrapper<ModelT, MemberFunctions::Objective>;

          std::fill_n(gradient, WrapperT::NUM_INPUTS, 0.0);

          __enzyme_autodiff<void>((void*) WrapperT::template eval<ScalarT>,
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
