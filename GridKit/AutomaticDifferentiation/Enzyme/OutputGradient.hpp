/**
 * @file OutputGradient.hpp
 * @brief Enzyme derivatives of computed model outputs.
 */
#pragma once

#include <array>

#include <GridKit/AutomaticDifferentiation/Enzyme/EnzymeDefinitions.hpp>

namespace GridKit
{
  namespace Enzyme
  {
    /// Differentiate a pure output function with respect to its local inputs.
    /// Computed outputs have no residual rows to sparsify; consumers compose
    /// these local derivatives through the existing signal Jacobian interface.
    template <typename ModelT, size_t N>
    struct OutputGradient
    {
      using ScalarT = typename ModelT::ScalarT;
      using Outputs = typename ModelT::Outputs;

      static ScalarT evaluate(const ModelT* model, Outputs output, const ScalarT* input)
      {
        return model->evaluateOutput(output, input);
      }

      static std::array<ScalarT, N> eval(const ModelT* model, Outputs output, const std::array<ScalarT, N>& input)
      {
        std::array<ScalarT, N> gradient{}, direction{};
        for (size_t n = 0; n < N; ++n)
        {
          direction.fill(ScalarT{0});
          direction[n] = ScalarT{1};
          gradient[n]  = Sparse::__enzyme_fwddiff<ScalarT>(
              (void*) evaluate, enzyme_const, model, enzyme_const, output, enzyme_dup, input.data(), direction.data());
        }
        return gradient;
      }
    };
  } // namespace Enzyme
} // namespace GridKit
