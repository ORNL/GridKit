#pragma once

#include <array>
#include <cstddef>

#include <GridKit/Optimization/Evaluator.hpp>
#include <GridKit/Optimization/SparseMatrixAssembler.hpp>

namespace GridKit
{
  namespace Testing
  {
    namespace Optimization
    {
      template <class scalar_type, typename index_type>
      class AffineModel : public GridKit::Optimization::Evaluator<scalar_type, index_type>
      {
      public:
        using ScalarT    = scalar_type;
        using IdxT       = index_type;
        using BaseT      = GridKit::Optimization::Evaluator<ScalarT, IdxT>;
        using RealT      = typename BaseT::RealT;
        using VectorT    = typename BaseT::VectorT;
        using CsrMatrixT = typename BaseT::CsrMatrixT;
        using AssemblerT = GridKit::Optimization::SparseMatrixAssembler<RealT, IdxT>;

        int allocate() override
        {
          allocated_ = false;

          int status  = variables_.resize(1);
          status     += variable_lower_bounds_.resize(1);
          status     += variable_upper_bounds_.resize(1);
          status     += objective_gradient_.resize(1);
          if (status != 0)
          {
            return status;
          }

          const std::array<typename AssemblerT::EntrySpan, 1> empty_groups{{typename AssemblerT::EntrySpan{}}};
          status      = jacobian_assembler_.allocate(0, 1, empty_groups, jacobian_contribution_);
          status     += hessian_assembler_.allocate(1, 1, empty_groups, hessian_contribution_);
          allocated_  = status == 0;
          return status;
        }

        int initialize() override
        {
          if (!allocated_)
          {
            return 1;
          }

          const std::array<ScalarT, 1> variables{{0.0}};
          const std::array<ScalarT, 1> lower_bounds{{-1.0}};
          const std::array<ScalarT, 1> upper_bounds{{1.0}};
          int                          status = variables_.copyFromExternal(
              variables.data(), memory::HOST, memory::HOST);
          status += variable_lower_bounds_.copyFromExternal(
              lower_bounds.data(), memory::HOST, memory::HOST);
          status += variable_upper_bounds_.copyFromExternal(
              upper_bounds.data(), memory::HOST, memory::HOST);
          status += objective_gradient_.setToZero(memory::HOST);
          return status;
        }

        IdxT size() override
        {
          return 1;
        }

        IdxT sizeConstraints() override
        {
          return 0;
        }

        VectorT& variables() override
        {
          return variables_;
        }

        const VectorT& variables() const override
        {
          return variables_;
        }

        const VectorT& variableLowerBounds() const override
        {
          return variable_lower_bounds_;
        }

        const VectorT& variableUpperBounds() const override
        {
          return variable_upper_bounds_;
        }

        int evaluateObjective() override
        {
          return 0;
        }

        RealT objective() const override
        {
          return RealT{0};
        }

        int evaluateObjectiveGradient() override
        {
          return objective_gradient_.setToZero(memory::HOST);
        }

        const VectorT& objectiveGradient() const override
        {
          return objective_gradient_;
        }

        int evaluateConstraints() override
        {
          return 0;
        }

        const VectorT& constraints() const override
        {
          return constraints_;
        }

        const VectorT& constraintLowerBounds() const override
        {
          return constraint_lower_bounds_;
        }

        const VectorT& constraintUpperBounds() const override
        {
          return constraint_upper_bounds_;
        }

        int evaluateJacobian() override
        {
          return jacobian_assembler_.clearValues();
        }

        CsrMatrixT* getCsrJacobian() const override
        {
          return jacobian_assembler_.matrix();
        }

        bool hasJacobian() override
        {
          return true;
        }

        int evaluateHessian(RealT, const RealT*, IdxT multiplier_count) override
        {
          if (multiplier_count != 0)
          {
            return 1;
          }
          const int status = hessian_assembler_.clearValues();
          if (status == 0)
          {
            ++hessian_evaluation_count_;
          }
          return status;
        }

        CsrMatrixT* getCsrHessian() const override
        {
          return hessian_assembler_.matrix();
        }

        bool hasHessian() override
        {
          return true;
        }

        std::size_t hessianEvaluationCount() const
        {
          return hessian_evaluation_count_;
        }

      private:
        VectorT variables_;
        VectorT variable_lower_bounds_;
        VectorT variable_upper_bounds_;
        VectorT objective_gradient_;
        VectorT constraints_;
        VectorT constraint_lower_bounds_;
        VectorT constraint_upper_bounds_;

        AssemblerT jacobian_assembler_;
        AssemblerT hessian_assembler_;

        std::array<typename AssemblerT::Contribution, 1> jacobian_contribution_;
        std::array<typename AssemblerT::Contribution, 1> hessian_contribution_;

        bool        allocated_{false};
        std::size_t hessian_evaluation_count_{0};
      };
    } // namespace Optimization
  } // namespace Testing
} // namespace GridKit
