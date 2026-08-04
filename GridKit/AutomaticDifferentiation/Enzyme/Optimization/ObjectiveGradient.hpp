#pragma once

#include <algorithm>
#include <cstddef>
#include <type_traits>

#include <GridKit/AutomaticDifferentiation/Enzyme/Optimization/ModelWrappers.hpp>

namespace GridKit
{
  namespace Enzyme
  {
    namespace Optimization
    {
      /**
       * @brief Exact reverse-mode objective gradient.
       */
      template <typename ModelT>
      struct ObjectiveGradient
      {
        using ScalarT = typename ModelT::ScalarT;

        static int eval(ModelT*        model,
                        std::size_t    variable_count,
                        const ScalarT* variables,
                        ScalarT*       gradient)
        {
          static_assert(std::is_same_v<ScalarT, double>,
                        "Optimization Enzyme derivatives currently support double only");

          if (model == nullptr
              || (variable_count > 0 && (variables == nullptr || gradient == nullptr)))
          {
            return 1;
          }

          std::fill_n(gradient, variable_count, ScalarT{0});

          GridKit::Enzyme::__enzyme_autodiff<void>(
              (void*) ModelWrapper<ModelT>::objective,
              enzyme_const,
              model,
              enzyme_dup,
              variables,
              gradient);
          return 0;
        }
      };
    } // namespace Optimization
  } // namespace Enzyme
} // namespace GridKit
