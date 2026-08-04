#pragma once

#include <algorithm>
#include <array>
#include <cstddef>

#include <GridKit/AutomaticDifferentiation/Enzyme/EnzymeDefinitions.hpp>

namespace GridKit
{
  namespace Enzyme
  {
    namespace Optimization
    {
      /**
       * @brief Primal model functions exposed to Enzyme.
       *
       * ModelT must provide a compile-time CONSTRAINT_COUNT and primal-only
       * evaluateObjective() and evaluateConstraints() methods. Their concrete,
       * nonvirtual definitions must be visible in the Enzyme-enabled
       * translation unit. Following the PhasorDynamics convention, include the
       * model's implementation header from its ModelEnzyme.cpp file. Enzyme
       * treats the model object as constant: model parameters may be read from
       * it, but every decision-variable-dependent load must use the variables
       * argument.
       */
      template <typename ModelT>
      struct ModelWrapper
      {
        using ScalarT = typename ModelT::ScalarT;

        __attribute__((always_inline)) inline static ScalarT
        objective(ModelT* model, const ScalarT* variables)
        {
          return model->evaluateObjective(variables);
        }

        __attribute__((always_inline)) inline static void
        constraints(ModelT* model, const ScalarT* variables, ScalarT* values)
        {
          model->evaluateConstraints(variables, values);
        }

        __attribute__((always_inline)) inline static ScalarT
        lagrangian(ModelT*        model,
                   const ScalarT* variables,
                   ScalarT        objective_factor,
                   const ScalarT* multipliers)
        {
          std::array<ScalarT, ModelT::CONSTRAINT_COUNT> constraint_values{};
          model->evaluateConstraints(variables, constraint_values.data());

          ScalarT value = objective_factor * model->evaluateObjective(variables);
          for (std::size_t row = 0; row < ModelT::CONSTRAINT_COUNT; ++row)
          {
            value += multipliers[row] * constraint_values[row];
          }
          return value;
        }

        __attribute__((always_inline)) inline static void
        lagrangianGradient(ModelT*        model,
                           std::size_t    variable_count,
                           const ScalarT* variables,
                           ScalarT        objective_factor,
                           const ScalarT* multipliers,
                           ScalarT*       gradient)
        {
          std::fill_n(gradient, variable_count, ScalarT{0});

          GridKit::Enzyme::__enzyme_autodiff<void>(
              (void*) lagrangian,
              enzyme_const,
              model,
              enzyme_dup,
              variables,
              gradient,
              enzyme_const,
              objective_factor,
              enzyme_const,
              multipliers);
        }
      };
    } // namespace Optimization
  } // namespace Enzyme
} // namespace GridKit
