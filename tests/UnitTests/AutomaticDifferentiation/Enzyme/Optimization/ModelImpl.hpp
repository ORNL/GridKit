#pragma once

#include <algorithm>
#include <array>
#include <vector>

#include "Model.hpp"

namespace GridKit
{
  namespace Testing
  {
    namespace Optimization
    {
      template <class scalar_type, typename index_type>
      void LeafA<scalar_type, index_type>::gatherVariables(const ScalarT* variables)
      {
        for (std::size_t entry = 0; entry < VARIABLE_COUNT; ++entry)
        {
          variables_[entry] = variables[static_cast<std::size_t>(VARIABLE_INDICES[entry])];
        }
      }

      template <class scalar_type, typename index_type>
      __attribute__((always_inline)) inline scalar_type
      LeafA<scalar_type, index_type>::evaluateObjective(const ScalarT* variables) const
      {
        const ScalarT x2 = variables[0];
        const ScalarT x0 = variables[1];
        return x2 * x2 * x2 / 3.0
               + x2 * x0
               + 0.5 * x0 * x0
               - 2.0 * x2
               - 2.0 * x0;
      }

      template <class scalar_type, typename index_type>
      __attribute__((always_inline)) inline void
      LeafA<scalar_type, index_type>::evaluateConstraints(
          const ScalarT* variables, ScalarT* constraints) const
      {
        const ScalarT x2 = variables[0];
        const ScalarT x0 = variables[1];
        constraints[0]   = x2 * x0 + 0.5 * x2 * x2;
      }

      template <class scalar_type, typename index_type>
      const std::array<scalar_type, LeafA<scalar_type, index_type>::VARIABLE_COUNT>&
      LeafA<scalar_type, index_type>::variables() const
      {
        return variables_;
      }

      template <class scalar_type, typename index_type>
      const std::array<scalar_type, LeafA<scalar_type, index_type>::VARIABLE_COUNT>&
      LeafA<scalar_type, index_type>::objectiveGradient() const
      {
        return objective_gradient_;
      }

      template <class scalar_type, typename index_type>
      const std::array<typename LeafA<scalar_type, index_type>::RealT,
                       LeafA<scalar_type, index_type>::JACOBIAN_ENTRIES.size()>&
      LeafA<scalar_type, index_type>::jacobianValues() const
      {
        return jacobian_values_;
      }

      template <class scalar_type, typename index_type>
      const std::array<typename LeafA<scalar_type, index_type>::RealT,
                       LeafA<scalar_type, index_type>::HESSIAN_ENTRIES.size()>&
      LeafA<scalar_type, index_type>::hessianValues() const
      {
        return hessian_values_;
      }

      template <class scalar_type, typename index_type>
      bool LeafA<scalar_type, index_type>::hasJacobian() const
      {
        return true;
      }

      template <class scalar_type, typename index_type>
      bool LeafA<scalar_type, index_type>::hasHessian() const
      {
        return true;
      }

      template <class scalar_type, typename index_type>
      void LeafB<scalar_type, index_type>::gatherVariables(const ScalarT* variables)
      {
        for (std::size_t entry = 0; entry < VARIABLE_COUNT; ++entry)
        {
          variables_[entry] = variables[static_cast<std::size_t>(VARIABLE_INDICES[entry])];
        }
      }

      template <class scalar_type, typename index_type>
      __attribute__((always_inline)) inline scalar_type
      LeafB<scalar_type, index_type>::evaluateObjective(const ScalarT* variables) const
      {
        const ScalarT x2 = variables[0];
        const ScalarT x1 = variables[1];
        return 0.5 * x2 * x2
               + x2 * x1
               + 0.5 * x1 * x1
               - 2.0 * x2
               - 2.0 * x1;
      }

      template <class scalar_type, typename index_type>
      __attribute__((always_inline)) inline void
      LeafB<scalar_type, index_type>::evaluateConstraints(
          const ScalarT* variables, ScalarT* constraints) const
      {
        const ScalarT x2 = variables[0];
        const ScalarT x1 = variables[1];
        constraints[0]   = x2 * x1 + 0.5 * x2 * x2;
        constraints[1]   = x2 + 0.5 * x1 * x1;
      }

      template <class scalar_type, typename index_type>
      const std::array<scalar_type, LeafB<scalar_type, index_type>::VARIABLE_COUNT>&
      LeafB<scalar_type, index_type>::variables() const
      {
        return variables_;
      }

      template <class scalar_type, typename index_type>
      const std::array<scalar_type, LeafB<scalar_type, index_type>::VARIABLE_COUNT>&
      LeafB<scalar_type, index_type>::objectiveGradient() const
      {
        return objective_gradient_;
      }

      template <class scalar_type, typename index_type>
      const std::array<typename LeafB<scalar_type, index_type>::RealT,
                       LeafB<scalar_type, index_type>::JACOBIAN_ENTRIES.size()>&
      LeafB<scalar_type, index_type>::jacobianValues() const
      {
        return jacobian_values_;
      }

      template <class scalar_type, typename index_type>
      const std::array<typename LeafB<scalar_type, index_type>::RealT,
                       LeafB<scalar_type, index_type>::HESSIAN_ENTRIES.size()>&
      LeafB<scalar_type, index_type>::hessianValues() const
      {
        return hessian_values_;
      }

      template <class scalar_type, typename index_type>
      bool LeafB<scalar_type, index_type>::hasJacobian() const
      {
        return true;
      }

      template <class scalar_type, typename index_type>
      bool LeafB<scalar_type, index_type>::hasHessian() const
      {
        return true;
      }

      template <class scalar_type, typename index_type>
      int Model<scalar_type, index_type>::allocate()
      {
        allocated_ = false;

        if (!GridKit::Optimization::hasUniqueIndices<IdxT>(
                LeafA<ScalarT, IdxT>::VARIABLE_INDICES)
            || !GridKit::Optimization::hasUniqueIndices<IdxT>(
                LeafB<ScalarT, IdxT>::VARIABLE_INDICES))
        {
          return 1;
        }

        int status  = 0;
        status     += variables_.resize(VARIABLE_COUNT);
        status     += variable_lower_bounds_.resize(VARIABLE_COUNT);
        status     += variable_upper_bounds_.resize(VARIABLE_COUNT);
        status     += objective_gradient_.resize(VARIABLE_COUNT);
        status     += constraints_.resize(CONSTRAINT_COUNT);
        status     += constraint_lower_bounds_.resize(CONSTRAINT_COUNT);
        status     += constraint_upper_bounds_.resize(CONSTRAINT_COUNT);
        if (status != 0)
        {
          return status;
        }

        std::vector<GridKit::Optimization::SparseEntry<IdxT>> jacobian_entries_a;
        std::vector<GridKit::Optimization::SparseEntry<IdxT>> jacobian_entries_b;
        jacobian_entries_a.reserve(LeafA<ScalarT, IdxT>::JACOBIAN_ENTRIES.size());
        jacobian_entries_b.reserve(LeafB<ScalarT, IdxT>::JACOBIAN_ENTRIES.size());

        for (const auto& entry : LeafA<ScalarT, IdxT>::JACOBIAN_ENTRIES)
        {
          jacobian_entries_a.push_back(
              {LeafA<ScalarT, IdxT>::CONSTRAINT_INDICES[static_cast<std::size_t>(entry.constraint)],
               LeafA<ScalarT, IdxT>::VARIABLE_INDICES[static_cast<std::size_t>(entry.variable)]});
        }
        for (const auto& entry : LeafB<ScalarT, IdxT>::JACOBIAN_ENTRIES)
        {
          jacobian_entries_b.push_back(
              {LeafB<ScalarT, IdxT>::CONSTRAINT_INDICES[static_cast<std::size_t>(entry.constraint)],
               LeafB<ScalarT, IdxT>::VARIABLE_INDICES[static_cast<std::size_t>(entry.variable)]});
        }

        const std::array<typename SparseAssemblerT::EntrySpan, 2> jacobian_groups{{jacobian_entries_a, jacobian_entries_b}};
        status = jacobian_assembler_.allocate(
            CONSTRAINT_COUNT,
            VARIABLE_COUNT,
            jacobian_groups,
            jacobian_contributions_);
        if (status != 0)
        {
          return status;
        }

        std::vector<GridKit::Optimization::SparseEntry<IdxT>> hessian_entries_a;
        std::vector<GridKit::Optimization::SparseEntry<IdxT>> hessian_entries_b;
        hessian_entries_a.reserve(LeafA<ScalarT, IdxT>::HESSIAN_ENTRIES.size());
        hessian_entries_b.reserve(LeafB<ScalarT, IdxT>::HESSIAN_ENTRIES.size());

        for (const auto& entry : LeafA<ScalarT, IdxT>::HESSIAN_ENTRIES)
        {
          hessian_entries_a.push_back(GridKit::Optimization::lowerTriangle(
              LeafA<ScalarT, IdxT>::VARIABLE_INDICES[static_cast<std::size_t>(entry.row)],
              LeafA<ScalarT, IdxT>::VARIABLE_INDICES[static_cast<std::size_t>(entry.column)]));
        }
        for (const auto& entry : LeafB<ScalarT, IdxT>::HESSIAN_ENTRIES)
        {
          hessian_entries_b.push_back(GridKit::Optimization::lowerTriangle(
              LeafB<ScalarT, IdxT>::VARIABLE_INDICES[static_cast<std::size_t>(entry.row)],
              LeafB<ScalarT, IdxT>::VARIABLE_INDICES[static_cast<std::size_t>(entry.column)]));
        }

        const std::array<typename SparseAssemblerT::EntrySpan, 2> hessian_groups{{hessian_entries_a, hessian_entries_b}};
        status = hessian_assembler_.allocate(
            VARIABLE_COUNT,
            VARIABLE_COUNT,
            hessian_groups,
            hessian_contributions_);
        if (status != 0)
        {
          return status;
        }

        allocated_ = true;
        return 0;
      }

      template <class scalar_type, typename index_type>
      int Model<scalar_type, index_type>::initialize()
      {
        if (!allocated_)
        {
          return 1;
        }

        auto* variable_values = variables_.getData(memory::HOST);
        variable_values[0]    = 0.8;
        variable_values[1]    = 1.2;
        variable_values[2]    = 0.9;
        variables_.setDataUpdated(memory::HOST);

        auto* lower_bounds = variable_lower_bounds_.getData(memory::HOST);
        auto* upper_bounds = variable_upper_bounds_.getData(memory::HOST);
        for (IdxT variable = 0; variable < VARIABLE_COUNT; ++variable)
        {
          lower_bounds[static_cast<std::size_t>(variable)] = 0.6;
          upper_bounds[static_cast<std::size_t>(variable)] = 1.4;
        }
        variable_lower_bounds_.setDataUpdated(memory::HOST);
        variable_upper_bounds_.setDataUpdated(memory::HOST);

        constraint_lower_bounds_.setToZero(memory::HOST);
        constraint_upper_bounds_.setToZero(memory::HOST);
        objective_gradient_.setToZero(memory::HOST);
        constraints_.setToZero(memory::HOST);
        return 0;
      }

      template <class scalar_type, typename index_type>
      index_type Model<scalar_type, index_type>::size()
      {
        return VARIABLE_COUNT;
      }

      template <class scalar_type, typename index_type>
      index_type Model<scalar_type, index_type>::sizeConstraints()
      {
        return CONSTRAINT_COUNT;
      }

      template <class scalar_type, typename index_type>
      typename Model<scalar_type, index_type>::VectorT&
      Model<scalar_type, index_type>::variables()
      {
        return variables_;
      }

      template <class scalar_type, typename index_type>
      const typename Model<scalar_type, index_type>::VectorT&
      Model<scalar_type, index_type>::variables() const
      {
        return variables_;
      }

      template <class scalar_type, typename index_type>
      const typename Model<scalar_type, index_type>::VectorT&
      Model<scalar_type, index_type>::variableLowerBounds() const
      {
        return variable_lower_bounds_;
      }

      template <class scalar_type, typename index_type>
      const typename Model<scalar_type, index_type>::VectorT&
      Model<scalar_type, index_type>::variableUpperBounds() const
      {
        return variable_upper_bounds_;
      }

      template <class scalar_type, typename index_type>
      int Model<scalar_type, index_type>::evaluateObjective()
      {
        if (gatherVariables() != 0)
        {
          return 1;
        }

        objective_ = leaf_a_.evaluateObjective(leaf_a_.variables().data())
                     + leaf_b_.evaluateObjective(leaf_b_.variables().data());
        return 0;
      }

      template <class scalar_type, typename index_type>
      typename Model<scalar_type, index_type>::RealT
      Model<scalar_type, index_type>::objective() const
      {
        return objective_;
      }

      template <class scalar_type, typename index_type>
      const typename Model<scalar_type, index_type>::VectorT&
      Model<scalar_type, index_type>::objectiveGradient() const
      {
        return objective_gradient_;
      }

      template <class scalar_type, typename index_type>
      int Model<scalar_type, index_type>::evaluateConstraints()
      {
        if (gatherVariables() != 0)
        {
          return 1;
        }

        std::array<ScalarT, LeafA<ScalarT, IdxT>::CONSTRAINT_COUNT> constraints_a{};
        std::array<ScalarT, LeafB<ScalarT, IdxT>::CONSTRAINT_COUNT> constraints_b{};
        leaf_a_.evaluateConstraints(leaf_a_.variables().data(), constraints_a.data());
        leaf_b_.evaluateConstraints(leaf_b_.variables().data(), constraints_b.data());

        auto* constraint_values = constraints_.getData(memory::HOST);
        constraint_values[0]    = constraints_a[0] + constraints_b[0] - 3.0;
        constraint_values[1]    = constraints_b[1] - 1.5;
        return constraints_.setDataUpdated(memory::HOST);
      }

      template <class scalar_type, typename index_type>
      const typename Model<scalar_type, index_type>::VectorT&
      Model<scalar_type, index_type>::constraints() const
      {
        return constraints_;
      }

      template <class scalar_type, typename index_type>
      const typename Model<scalar_type, index_type>::VectorT&
      Model<scalar_type, index_type>::constraintLowerBounds() const
      {
        return constraint_lower_bounds_;
      }

      template <class scalar_type, typename index_type>
      const typename Model<scalar_type, index_type>::VectorT&
      Model<scalar_type, index_type>::constraintUpperBounds() const
      {
        return constraint_upper_bounds_;
      }

      template <class scalar_type, typename index_type>
      typename Model<scalar_type, index_type>::CsrMatrixT*
      Model<scalar_type, index_type>::getCsrJacobian() const
      {
        return jacobian_assembler_.matrix();
      }

      template <class scalar_type, typename index_type>
      typename Model<scalar_type, index_type>::CsrMatrixT*
      Model<scalar_type, index_type>::getCsrHessian() const
      {
        return hessian_assembler_.matrix();
      }

      template <class scalar_type, typename index_type>
      const LinearLeaf<index_type>& Model<scalar_type, index_type>::linearLeaf() const
      {
        return linear_leaf_;
      }

      template <class scalar_type, typename index_type>
      std::size_t Model<scalar_type, index_type>::hessianEvaluationCount() const
      {
        return hessian_evaluation_count_;
      }

      template <class scalar_type, typename index_type>
      int Model<scalar_type, index_type>::gatherVariables()
      {
        if (!allocated_)
        {
          return 1;
        }

        const auto* variables = variables_.getData(memory::HOST);
        leaf_a_.gatherVariables(variables);
        leaf_b_.gatherVariables(variables);
        return 0;
      }
    } // namespace Optimization
  } // namespace Testing
} // namespace GridKit
