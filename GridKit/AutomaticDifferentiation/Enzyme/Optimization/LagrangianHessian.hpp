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
       * @brief Exact fixed-pattern lower-triangular Lagrangian Hessian values.
       *
       * Entries must be lower triangular and ordered by column and then row.
       * Forward-over-reverse evaluates one Hessian-vector product for each
       * structural column; only requested entries are extracted.
       */
      template <typename ModelT, typename index_type>
      struct LagrangianHessian
      {
        using ScalarT = typename ModelT::ScalarT;
        using RealT   = typename ModelT::RealT;
        using IdxT    = index_type;
        using EntryT  = GridKit::Optimization::LocalHessianEntry<IdxT>;

        static int eval(ModelT*                 model,
                        std::size_t             variable_count,
                        std::span<const EntryT> entries,
                        const ScalarT*          variables,
                        ScalarT                 objective_factor,
                        const ScalarT*          multipliers,
                        RealT*                  values)
        {
          static_assert(std::is_same_v<ScalarT, double>,
                        "Optimization Enzyme derivatives currently support double only");

          if (entries.empty())
          {
            return 0;
          }
          if (model == nullptr || variables == nullptr || values == nullptr)
          {
            return 1;
          }
          if constexpr (ModelT::CONSTRAINT_COUNT > 0)
          {
            if (multipliers == nullptr)
            {
              return 1;
            }
          }
          for (std::size_t entry = 0; entry < entries.size(); ++entry)
          {
            const auto column = static_cast<std::size_t>(entries[entry].column);
            const auto row    = static_cast<std::size_t>(entries[entry].row);
            if (column >= variable_count || row >= variable_count || row < column)
            {
              return 1;
            }
            if (entry > 0)
            {
              const auto previous_column =
                  static_cast<std::size_t>(entries[entry - 1].column);
              const auto previous_row =
                  static_cast<std::size_t>(entries[entry - 1].row);
              if (column < previous_column
                  || (column == previous_column && row <= previous_row))
              {
                return 1;
              }
            }
          }

          std::vector<ScalarT> seed(variable_count);
          std::vector<ScalarT> primal_gradient(variable_count);
          std::vector<ScalarT> hessian_product(variable_count);

          std::size_t entry = 0;
          while (entry < entries.size())
          {
            const auto column = static_cast<std::size_t>(entries[entry].column);
            std::fill(seed.begin(), seed.end(), ScalarT{0});
            std::fill(primal_gradient.begin(), primal_gradient.end(), ScalarT{0});
            std::fill(hessian_product.begin(), hessian_product.end(), ScalarT{0});
            seed[column] = ScalarT{1};

            GridKit::Enzyme::__enzyme_fwddiff<void>(
                (void*) ModelWrapper<ModelT>::lagrangianGradient,
                enzyme_const,
                model,
                enzyme_const,
                variable_count,
                enzyme_dup,
                variables,
                seed.data(),
                enzyme_const,
                objective_factor,
                enzyme_const,
                multipliers,
                enzyme_dupnoneed,
                primal_gradient.data(),
                hessian_product.data());

            do
            {
              const auto row = static_cast<std::size_t>(entries[entry].row);
              values[entry]  = static_cast<RealT>(hessian_product[row]);
              ++entry;
            } while (entry < entries.size()
                     && static_cast<std::size_t>(entries[entry].column) == column);
          }
          return 0;
        }
      };
    } // namespace Optimization
  } // namespace Enzyme
} // namespace GridKit
