#pragma once

#include <array>
#include <cstddef>
#include <limits>
#include <type_traits>

#include <GridKit/Testing/Testing.hpp>

#include "Model.hpp"

namespace GridKit
{
  namespace Testing
  {
    namespace Optimization
    {
      template <class scalar_type, typename index_type>
      class OptimizationTests
      {
      public:
        using ScalarT    = scalar_type;
        using IdxT       = index_type;
        using ModelT     = Model<ScalarT, IdxT>;
        using RealT      = typename ModelT::RealT;
        using CsrMatrixT = typename ModelT::CsrMatrixT;

        TestOutcome exactSparseDerivatives()
        {
          TestStatus success = true;
          ModelT     model;

          success *= model.allocate() == 0;
          success *= model.initialize() == 0;
          success *= model.hasJacobian();
          success *= model.hasHessian();
          success *= model.linearLeaf().hasHessian();
          success *= model.linearLeaf().hessianEntries().empty();

          auto* jacobian  = model.getCsrJacobian();
          auto* hessian   = model.getCsrHessian();
          success        *= jacobian != nullptr;
          success        *= hessian != nullptr;
          if (!success)
          {
            return success.report(__func__);
          }

          const auto* jacobian_rows    = jacobian->getRowData(memory::HOST);
          const auto* jacobian_columns = jacobian->getColData(memory::HOST);
          const auto* hessian_rows     = hessian->getRowData(memory::HOST);
          const auto* hessian_columns  = hessian->getColData(memory::HOST);

          const std::array<ScalarT, 3> variables{{2.0, 3.0, 4.0}};
          success *= setVariables(model, variables);
          success *= model.evaluateObjective() == 0;
          success *= isEqual(model.objective(), static_cast<RealT>(179.0 / 6.0), tolerance());

          success *= model.evaluateObjectiveGradient() == 0;
          success *= checkVector(model.objectiveGradient(),
                                 std::array<RealT, 3>{{4.0, 5.0, 21.0}});

          success *= model.evaluateConstraints() == 0;
          success *= checkVector(model.constraints(),
                                 std::array<RealT, 2>{{33.0, 7.0}});

          success *= model.evaluateJacobian() == 0;
          success *= checkCsr(jacobian,
                              std::array<IdxT, 3>{{0, 3, 5}},
                              std::array<IdxT, 5>{{0, 1, 2, 1, 2}},
                              std::array<RealT, 5>{{4.0, 4.0, 13.0, 3.0, 1.0}});

          success *= checkHessian(model, 1.0, {0.0, 0.0}, {1.0, 1.0, 1.0, 1.0, 9.0});
          success *= checkHessian(model, 0.0, {1.0, 0.0}, {0.0, 0.0, 1.0, 1.0, 2.0});
          success *= checkHessian(model, 0.0, {0.0, 1.0}, {0.0, 1.0, 0.0, 0.0, 0.0});
          success *= checkHessian(model, 2.0, {5.0, 7.0}, {2.0, 9.0, 7.0, 7.0, 28.0});

          const std::array<RealT, 2> multipliers{{0.0, 0.0}};
          success *= model.evaluateHessian(1.0, multipliers.data(), 1) != 0;
          success *= model.evaluateHessian(1.0, nullptr, ModelT::CONSTRAINT_COUNT) != 0;

          const std::array<ScalarT, 3> zero_variables{{2.0, 0.0, 0.0}};
          const std::array<RealT, 2>   zero_multipliers{{-1.0, -2.0}};
          success *= setVariables(model, zero_variables);
          success *= model.evaluateJacobian() == 0;
          success *= model.evaluateHessian(2.0, zero_multipliers.data(), zero_multipliers.size()) == 0;

          success *= checkCsr(jacobian,
                              std::array<IdxT, 3>{{0, 3, 5}},
                              std::array<IdxT, 5>{{0, 1, 2, 1, 2}},
                              std::array<RealT, 5>{{0.0, 0.0, 2.0, 0.0, 1.0}});
          success *= checkCsr(hessian,
                              std::array<IdxT, 4>{{0, 1, 2, 5}},
                              std::array<IdxT, 5>{{0, 1, 0, 1, 2}},
                              std::array<RealT, 5>{{2.0, 0.0, 1.0, 1.0, 0.0}});

          success *= jacobian_rows == jacobian->getRowData(memory::HOST);
          success *= jacobian_columns == jacobian->getColData(memory::HOST);
          success *= hessian_rows == hessian->getRowData(memory::HOST);
          success *= hessian_columns == hessian->getColData(memory::HOST);

          const std::array<IdxT, 2> duplicate_map{{1, 1}};
          success *= !GridKit::Optimization::hasUniqueIndices<IdxT>(duplicate_map);

          using AssemblerT = GridKit::Optimization::SparseMatrixAssembler<RealT, IdxT>;
          static_assert(!std::is_copy_constructible_v<AssemblerT>);
          static_assert(!std::is_move_constructible_v<AssemblerT>);

          AssemblerT                                          empty_assembler;
          const std::array<typename AssemblerT::EntrySpan, 1> empty_groups{{typename AssemblerT::EntrySpan{}}};
          std::array<typename AssemblerT::Contribution, 1>    empty_contributions;
          success *= empty_assembler.allocate(3, 3, empty_groups, empty_contributions) == 0;
          success *= empty_assembler.matrix() != nullptr;
          if (empty_assembler.matrix() != nullptr)
          {
            success          *= empty_assembler.matrix()->getNnz() == 0;
            success          *= empty_contributions[0].size() == 0;
            success          *= empty_assembler.clearValues() == 0;
            success          *= empty_assembler.addValues(empty_contributions[0], {}) == 0;
            const auto* rows  = empty_assembler.matrix()->getRowData(memory::HOST);
            success          *= rows[0] == 0 && rows[1] == 0 && rows[2] == 0 && rows[3] == 0;

            const auto stale_contribution  = empty_contributions[0];
            success                       *= empty_assembler.allocate(
                           3, 3, empty_groups, empty_contributions)
                       == 0;
            success *= empty_assembler.addValues(stale_contribution, {}) != 0;

            const auto* matrix_before_failure = empty_assembler.matrix();
            const std::array<GridKit::Optimization::SparseEntry<IdxT>, 1>
                                                                invalid_entries{{{3, 0}}};
            const std::array<typename AssemblerT::EntrySpan, 1> invalid_groups{{invalid_entries}};
            std::array<typename AssemblerT::Contribution, 1>    invalid_contributions;
            success *= empty_assembler.allocate(
                           3, 3, invalid_groups, invalid_contributions)
                       != 0;
            success *= empty_assembler.matrix() == matrix_before_failure;
            success *= empty_assembler.addValues(empty_contributions[0], {}) == 0;
          }

          return success.report(__func__);
        }

      private:
        static constexpr RealT tolerance()
        {
          return static_cast<RealT>(100.0) * std::numeric_limits<RealT>::epsilon();
        }

        template <std::size_t Count>
        static bool setVariables(ModelT& model, const std::array<ScalarT, Count>& values)
        {
          if (Count != static_cast<std::size_t>(model.size()))
          {
            return false;
          }
          return model.variables().copyFromExternal(
                     values.data(), memory::HOST, memory::HOST)
                 == 0;
        }

        template <std::size_t Count>
        static bool checkVector(const typename ModelT::VectorT& vector,
                                const std::array<RealT, Count>& expected,
                                RealT                           comparison_tolerance = tolerance())
        {
          if (vector.getSize() != static_cast<IdxT>(Count))
          {
            return false;
          }

          const auto* values = vector.getData(memory::HOST);
          for (std::size_t entry = 0; entry < Count; ++entry)
          {
            if (!isEqual(static_cast<RealT>(values[entry]),
                         expected[entry],
                         comparison_tolerance))
            {
              return false;
            }
          }
          return true;
        }

        template <std::size_t RowCount, std::size_t NonzeroCount>
        static bool checkCsr(CsrMatrixT*                            matrix,
                             const std::array<IdxT, RowCount>&      expected_rows,
                             const std::array<IdxT, NonzeroCount>&  expected_columns,
                             const std::array<RealT, NonzeroCount>& expected_values)
        {
          if (matrix == nullptr
              || matrix->getNumRows() + 1 != static_cast<IdxT>(RowCount)
              || matrix->getNumColumns() != ModelT::VARIABLE_COUNT
              || matrix->getNnz() != static_cast<IdxT>(NonzeroCount))
          {
            return false;
          }

          const auto* rows    = matrix->getRowData(memory::HOST);
          const auto* columns = matrix->getColData(memory::HOST);
          const auto* values  = matrix->getValues(memory::HOST);

          for (std::size_t entry = 0; entry < RowCount; ++entry)
          {
            if (rows[entry] != expected_rows[entry])
            {
              return false;
            }
          }
          for (std::size_t entry = 0; entry < NonzeroCount; ++entry)
          {
            if (columns[entry] != expected_columns[entry]
                || !isEqual(values[entry], expected_values[entry], tolerance()))
            {
              return false;
            }
          }
          return true;
        }

        static bool checkHessian(ModelT&              model,
                                 RealT                objective_factor,
                                 std::array<RealT, 2> multipliers,
                                 std::array<RealT, 5> expected_values)
        {
          if (model.evaluateHessian(
                  objective_factor, multipliers.data(), multipliers.size())
              != 0)
          {
            return false;
          }
          return checkCsr(model.getCsrHessian(),
                          std::array<IdxT, 4>{{0, 1, 2, 5}},
                          std::array<IdxT, 5>{{0, 1, 0, 1, 2}},
                          expected_values);
        }
      };
    } // namespace Optimization
  } // namespace Testing
} // namespace GridKit
