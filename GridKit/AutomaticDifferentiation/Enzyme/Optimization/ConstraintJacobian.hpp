#pragma once

#include <algorithm>
#include <cstddef>
#include <span>
#include <type_traits>
#include <vector>

#include <GridKit/AutomaticDifferentiation/Enzyme/Optimization/ModelWrappers.hpp>
#include <GridKit/Optimization/DerivativeStructure.hpp>

namespace GridKit
{
  namespace Enzyme
  {
    namespace Optimization
    {
      /**
       * @brief Exact fixed-pattern constraint Jacobian values.
       *
       * Entries must be ordered by variable and then constraint. One exact
       * Jacobian-vector product is evaluated for each structural column.
       */
      template <typename ModelT, typename index_type>
      struct ConstraintJacobian
      {
        using ScalarT = typename ModelT::ScalarT;
        using RealT   = typename ModelT::RealT;
        using IdxT    = index_type;
        using EntryT  = GridKit::Optimization::LocalJacobianEntry<IdxT>;

        static int eval(ModelT*                 model,
                        std::size_t             variable_count,
                        std::span<const EntryT> entries,
                        const ScalarT*          variables,
                        RealT*                  values)
        {
          static_assert(std::is_same_v<ScalarT, double>,
                        "Optimization Enzyme derivatives currently support double only");

          constexpr std::size_t constraint_count = ModelT::CONSTRAINT_COUNT;

          if (entries.empty())
          {
            return 0;
          }
          if (model == nullptr || variables == nullptr || values == nullptr)
          {
            return 1;
          }
          for (std::size_t entry = 0; entry < entries.size(); ++entry)
          {
            const auto column = static_cast<std::size_t>(entries[entry].variable);
            const auto row    = static_cast<std::size_t>(entries[entry].constraint);
            if (column >= variable_count || row >= constraint_count)
            {
              return 1;
            }
            if (entry > 0)
            {
              const auto previous_column =
                  static_cast<std::size_t>(entries[entry - 1].variable);
              const auto previous_row =
                  static_cast<std::size_t>(entries[entry - 1].constraint);
              if (column < previous_column
                  || (column == previous_column && row <= previous_row))
              {
                return 1;
              }
            }
          }

          std::vector<ScalarT> seed(variable_count);
          std::vector<ScalarT> primal_constraints(constraint_count);
          std::vector<ScalarT> constraint_tangent(constraint_count);

          std::size_t entry = 0;
          while (entry < entries.size())
          {
            const auto column = static_cast<std::size_t>(entries[entry].variable);
            std::fill(seed.begin(), seed.end(), ScalarT{0});
            std::fill(primal_constraints.begin(), primal_constraints.end(), ScalarT{0});
            std::fill(constraint_tangent.begin(), constraint_tangent.end(), ScalarT{0});
            seed[column] = ScalarT{1};

            GridKit::Enzyme::__enzyme_fwddiff<void>(
                (void*) ModelWrapper<ModelT>::constraints,
                enzyme_const,
                model,
                enzyme_dup,
                variables,
                seed.data(),
                enzyme_dupnoneed,
                primal_constraints.data(),
                constraint_tangent.data());

            do
            {
              const auto row = static_cast<std::size_t>(entries[entry].constraint);
              values[entry]  = static_cast<RealT>(constraint_tangent[row]);
              ++entry;
            } while (entry < entries.size()
                     && static_cast<std::size_t>(entries[entry].variable) == column);
          }
          return 0;
        }
      };
    } // namespace Optimization
  } // namespace Enzyme
} // namespace GridKit
