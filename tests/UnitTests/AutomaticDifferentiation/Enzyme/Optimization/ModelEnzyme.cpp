#include <GridKit/AutomaticDifferentiation/Enzyme/Optimization/SparseDerivatives.hpp>

#include "ModelImpl.hpp"

namespace GridKit
{
  namespace Testing
  {
    namespace Optimization
    {
      template <class scalar_type, typename index_type>
      int LeafA<scalar_type, index_type>::evaluateObjectiveGradient()
      {
        return GridKit::Enzyme::Optimization::ObjectiveGradient<LeafA>::eval(
            this, VARIABLE_COUNT, variables_.data(), objective_gradient_.data());
      }

      template <class scalar_type, typename index_type>
      int LeafA<scalar_type, index_type>::evaluateJacobian()
      {
        return GridKit::Enzyme::Optimization::ConstraintJacobian<LeafA, IdxT>::eval(
            this,
            VARIABLE_COUNT,
            JACOBIAN_ENTRIES,
            variables_.data(),
            jacobian_values_.data());
      }

      template <class scalar_type, typename index_type>
      int LeafA<scalar_type, index_type>::evaluateHessian(
          RealT objective_factor, const RealT* multipliers)
      {
        return GridKit::Enzyme::Optimization::LagrangianHessian<LeafA, IdxT>::eval(
            this,
            VARIABLE_COUNT,
            HESSIAN_ENTRIES,
            variables_.data(),
            objective_factor,
            multipliers,
            hessian_values_.data());
      }

      template <class scalar_type, typename index_type>
      int LeafB<scalar_type, index_type>::evaluateObjectiveGradient()
      {
        return GridKit::Enzyme::Optimization::ObjectiveGradient<LeafB>::eval(
            this, VARIABLE_COUNT, variables_.data(), objective_gradient_.data());
      }

      template <class scalar_type, typename index_type>
      int LeafB<scalar_type, index_type>::evaluateJacobian()
      {
        return GridKit::Enzyme::Optimization::ConstraintJacobian<LeafB, IdxT>::eval(
            this,
            VARIABLE_COUNT,
            JACOBIAN_ENTRIES,
            variables_.data(),
            jacobian_values_.data());
      }

      template <class scalar_type, typename index_type>
      int LeafB<scalar_type, index_type>::evaluateHessian(
          RealT objective_factor, const RealT* multipliers)
      {
        return GridKit::Enzyme::Optimization::LagrangianHessian<LeafB, IdxT>::eval(
            this,
            VARIABLE_COUNT,
            HESSIAN_ENTRIES,
            variables_.data(),
            objective_factor,
            multipliers,
            hessian_values_.data());
      }

      template <class scalar_type, typename index_type>
      int Model<scalar_type, index_type>::evaluateObjectiveGradient()
      {
        if (gatherVariables() != 0
            || leaf_a_.evaluateObjectiveGradient() != 0
            || leaf_b_.evaluateObjectiveGradient() != 0)
        {
          return 1;
        }

        objective_gradient_.setToZero(memory::HOST);
        auto* gradient = objective_gradient_.getData(memory::HOST);

        for (std::size_t entry = 0; entry < LeafA<ScalarT, IdxT>::VARIABLE_COUNT; ++entry)
        {
          const auto global                           = LeafA<ScalarT, IdxT>::VARIABLE_INDICES[entry];
          gradient[static_cast<std::size_t>(global)] += leaf_a_.objectiveGradient()[entry];
        }
        for (std::size_t entry = 0; entry < LeafB<ScalarT, IdxT>::VARIABLE_COUNT; ++entry)
        {
          const auto global                           = LeafB<ScalarT, IdxT>::VARIABLE_INDICES[entry];
          gradient[static_cast<std::size_t>(global)] += leaf_b_.objectiveGradient()[entry];
        }
        return objective_gradient_.setDataUpdated(memory::HOST);
      }

      template <class scalar_type, typename index_type>
      int Model<scalar_type, index_type>::evaluateJacobian()
      {
        if (gatherVariables() != 0
            || leaf_a_.evaluateJacobian() != 0
            || leaf_b_.evaluateJacobian() != 0
            || jacobian_assembler_.clearValues() != 0)
        {
          return 1;
        }

        if (jacobian_assembler_.addValues(
                jacobian_contributions_[0], leaf_a_.jacobianValues())
            != 0)
        {
          return 1;
        }
        return jacobian_assembler_.addValues(
            jacobian_contributions_[1], leaf_b_.jacobianValues());
      }

      template <class scalar_type, typename index_type>
      bool Model<scalar_type, index_type>::hasJacobian()
      {
        return leaf_a_.hasJacobian()
               && leaf_b_.hasJacobian()
               && linear_leaf_.hasJacobian();
      }

      template <class scalar_type, typename index_type>
      int Model<scalar_type, index_type>::evaluateHessian(
          RealT        objective_factor,
          const RealT* multipliers,
          IdxT         multiplier_count)
      {
        if (multiplier_count != CONSTRAINT_COUNT || multipliers == nullptr)
        {
          return 1;
        }
        if (gatherVariables() != 0)
        {
          return 1;
        }

        const std::array<RealT, LeafA<ScalarT, IdxT>::CONSTRAINT_COUNT> multipliers_a{{multipliers[0]}};
        const std::array<RealT, LeafB<ScalarT, IdxT>::CONSTRAINT_COUNT> multipliers_b{{multipliers[0], multipliers[1]}};

        if (leaf_a_.evaluateHessian(objective_factor, multipliers_a.data()) != 0
            || leaf_b_.evaluateHessian(objective_factor, multipliers_b.data()) != 0
            || hessian_assembler_.clearValues() != 0)
        {
          return 1;
        }

        if (hessian_assembler_.addValues(
                hessian_contributions_[0], leaf_a_.hessianValues())
            != 0)
        {
          return 1;
        }
        const int status = hessian_assembler_.addValues(
            hessian_contributions_[1], leaf_b_.hessianValues());
        if (status == 0)
        {
          ++hessian_evaluation_count_;
        }
        return status;
      }

      template <class scalar_type, typename index_type>
      bool Model<scalar_type, index_type>::hasHessian()
      {
        return leaf_a_.hasHessian()
               && leaf_b_.hasHessian()
               && linear_leaf_.hasHessian();
      }

      template class LeafA<double, std::size_t>;
      template class LeafB<double, std::size_t>;
      template class Model<double, std::size_t>;
    } // namespace Optimization
  } // namespace Testing
} // namespace GridKit
