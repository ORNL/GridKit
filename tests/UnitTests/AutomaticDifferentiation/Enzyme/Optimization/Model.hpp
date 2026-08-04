#pragma once

#include <array>
#include <cstddef>
#include <span>
#include <vector>

#include <GridKit/Optimization/DerivativeStructure.hpp>
#include <GridKit/Optimization/Evaluator.hpp>
#include <GridKit/Optimization/SparseMatrixAssembler.hpp>

namespace GridKit
{
  namespace Testing
  {
    namespace Optimization
    {
      template <class scalar_type, typename index_type>
      class LeafA
      {
      public:
        using ScalarT = scalar_type;
        using IdxT    = index_type;
        using RealT   = typename ScalarTraits<ScalarT>::RealT;

        static constexpr std::size_t VARIABLE_COUNT   = 2;
        static constexpr std::size_t CONSTRAINT_COUNT = 1;

        inline static constexpr std::array<IdxT, VARIABLE_COUNT>   VARIABLE_INDICES{{2, 0}};
        inline static constexpr std::array<IdxT, CONSTRAINT_COUNT> CONSTRAINT_INDICES{{0}};

        inline static constexpr std::array<GridKit::Optimization::LocalJacobianEntry<IdxT>, 2>
            JACOBIAN_ENTRIES{{{0, 0}, {0, 1}}};
        inline static constexpr std::array<GridKit::Optimization::LocalHessianEntry<IdxT>, 3>
            HESSIAN_ENTRIES{{{0, 0}, {1, 0}, {1, 1}}};

        void gatherVariables(const ScalarT* variables);

        __attribute__((always_inline)) inline ScalarT
        evaluateObjective(const ScalarT* variables) const;

        __attribute__((always_inline)) inline void
        evaluateConstraints(const ScalarT* variables, ScalarT* constraints) const;

        int evaluateObjectiveGradient();
        int evaluateJacobian();
        int evaluateHessian(RealT objective_factor, const RealT* multipliers);

        const std::array<ScalarT, VARIABLE_COUNT>&        variables() const;
        const std::array<ScalarT, VARIABLE_COUNT>&        objectiveGradient() const;
        const std::array<RealT, JACOBIAN_ENTRIES.size()>& jacobianValues() const;
        const std::array<RealT, HESSIAN_ENTRIES.size()>&  hessianValues() const;

        bool hasJacobian() const;
        bool hasHessian() const;

      private:
        std::array<ScalarT, VARIABLE_COUNT>        variables_{};
        std::array<ScalarT, VARIABLE_COUNT>        objective_gradient_{};
        std::array<RealT, JACOBIAN_ENTRIES.size()> jacobian_values_{};
        std::array<RealT, HESSIAN_ENTRIES.size()>  hessian_values_{};
      };

      template <class scalar_type, typename index_type>
      class LeafB
      {
      public:
        using ScalarT = scalar_type;
        using IdxT    = index_type;
        using RealT   = typename ScalarTraits<ScalarT>::RealT;

        static constexpr std::size_t VARIABLE_COUNT   = 2;
        static constexpr std::size_t CONSTRAINT_COUNT = 2;

        inline static constexpr std::array<IdxT, VARIABLE_COUNT>   VARIABLE_INDICES{{2, 1}};
        inline static constexpr std::array<IdxT, CONSTRAINT_COUNT> CONSTRAINT_INDICES{{0, 1}};

        inline static constexpr std::array<GridKit::Optimization::LocalJacobianEntry<IdxT>, 4>
            JACOBIAN_ENTRIES{{{0, 0}, {1, 0}, {0, 1}, {1, 1}}};
        inline static constexpr std::array<GridKit::Optimization::LocalHessianEntry<IdxT>, 3>
            HESSIAN_ENTRIES{{{0, 0}, {1, 0}, {1, 1}}};

        void gatherVariables(const ScalarT* variables);

        __attribute__((always_inline)) inline ScalarT
        evaluateObjective(const ScalarT* variables) const;

        __attribute__((always_inline)) inline void
        evaluateConstraints(const ScalarT* variables, ScalarT* constraints) const;

        int evaluateObjectiveGradient();
        int evaluateJacobian();
        int evaluateHessian(RealT objective_factor, const RealT* multipliers);

        const std::array<ScalarT, VARIABLE_COUNT>&        variables() const;
        const std::array<ScalarT, VARIABLE_COUNT>&        objectiveGradient() const;
        const std::array<RealT, JACOBIAN_ENTRIES.size()>& jacobianValues() const;
        const std::array<RealT, HESSIAN_ENTRIES.size()>&  hessianValues() const;

        bool hasJacobian() const;
        bool hasHessian() const;

      private:
        std::array<ScalarT, VARIABLE_COUNT>        variables_{};
        std::array<ScalarT, VARIABLE_COUNT>        objective_gradient_{};
        std::array<RealT, JACOBIAN_ENTRIES.size()> jacobian_values_{};
        std::array<RealT, HESSIAN_ENTRIES.size()>  hessian_values_{};
      };

      template <typename index_type>
      class LinearLeaf
      {
      public:
        using IdxT = index_type;

        std::span<const GridKit::Optimization::LocalHessianEntry<IdxT>> hessianEntries() const
        {
          return HESSIAN_ENTRIES;
        }

        bool hasJacobian() const
        {
          return true;
        }

        bool hasHessian() const
        {
          return true;
        }

      private:
        inline static constexpr std::array<GridKit::Optimization::LocalHessianEntry<IdxT>, 0>
            HESSIAN_ENTRIES{};
      };

      template <class scalar_type, typename index_type>
      class Model : public GridKit::Optimization::Evaluator<scalar_type, index_type>
      {
      public:
        using ScalarT    = scalar_type;
        using IdxT       = index_type;
        using BaseT      = GridKit::Optimization::Evaluator<ScalarT, IdxT>;
        using RealT      = typename BaseT::RealT;
        using VectorT    = typename BaseT::VectorT;
        using CsrMatrixT = typename BaseT::CsrMatrixT;
        using SparseAssemblerT =
            GridKit::Optimization::SparseMatrixAssembler<RealT, IdxT>;
        using ContributionT = typename SparseAssemblerT::Contribution;

        static constexpr IdxT VARIABLE_COUNT   = 3;
        static constexpr IdxT CONSTRAINT_COUNT = 2;

        int allocate() override;
        int initialize() override;

        IdxT size() override;
        IdxT sizeConstraints() override;

        VectorT&       variables() override;
        const VectorT& variables() const override;
        const VectorT& variableLowerBounds() const override;
        const VectorT& variableUpperBounds() const override;

        int   evaluateObjective() override;
        RealT objective() const override;

        int            evaluateObjectiveGradient() override;
        const VectorT& objectiveGradient() const override;

        int            evaluateConstraints() override;
        const VectorT& constraints() const override;
        const VectorT& constraintLowerBounds() const override;
        const VectorT& constraintUpperBounds() const override;

        int         evaluateJacobian() override;
        CsrMatrixT* getCsrJacobian() const override;
        bool        hasJacobian() override;

        int         evaluateHessian(RealT        objective_factor,
                                    const RealT* multipliers,
                                    IdxT         multiplier_count) override;
        CsrMatrixT* getCsrHessian() const override;
        bool        hasHessian() override;

        const LinearLeaf<IdxT>& linearLeaf() const;
        std::size_t             hessianEvaluationCount() const;

      private:
        int gatherVariables();

        VectorT variables_;
        VectorT variable_lower_bounds_;
        VectorT variable_upper_bounds_;
        VectorT objective_gradient_;
        VectorT constraints_;
        VectorT constraint_lower_bounds_;
        VectorT constraint_upper_bounds_;

        RealT objective_{0};

        LeafA<ScalarT, IdxT> leaf_a_;
        LeafB<ScalarT, IdxT> leaf_b_;
        LinearLeaf<IdxT>     linear_leaf_;

        SparseAssemblerT jacobian_assembler_;
        SparseAssemblerT hessian_assembler_;

        std::array<ContributionT, 2> jacobian_contributions_;
        std::array<ContributionT, 2> hessian_contributions_;

        bool        allocated_{false};
        std::size_t hessian_evaluation_count_{0};
      };

      extern template class Model<double, std::size_t>;
    } // namespace Optimization
  } // namespace Testing
} // namespace GridKit
