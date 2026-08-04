#pragma once

#include <cstddef>
#include <memory>
#include <vector>

#include <GridKit/Model/OPF/Component.hpp>
#include <GridKit/Model/OPF/SystemData.hpp>
#include <GridKit/Model/StateData.hpp>
#include <GridKit/Optimization/Evaluator.hpp>
#include <GridKit/Optimization/SparseMatrixAssembler.hpp>

namespace GridKit
{
  namespace OPF
  {
    /**
     * @brief Polar AC optimal-power-flow evaluator with exact sparse derivatives.
     *
     * Component topology and derivative structures are fixed by allocate().
     * Objective gradients, constraint Jacobians, and Lagrangian Hessians are
     * evaluated locally and scatter-added into system-owned storage.
     */
    template <class scalar_type = double, typename index_type = std::size_t>
    class System final : public Optimization::Evaluator<scalar_type, index_type>
    {
    public:
      using ScalarT       = scalar_type;
      using IdxT          = index_type;
      using BaseT         = Optimization::Evaluator<ScalarT, IdxT>;
      using RealT         = typename BaseT::RealT;
      using VectorT       = typename BaseT::VectorT;
      using CsrMatrixT    = typename BaseT::CsrMatrixT;
      using ComponentT    = Component<ScalarT, IdxT>;
      using SystemDataT   = SystemData<RealT, IdxT>;
      using AssemblerT    = Optimization::SparseMatrixAssembler<RealT, IdxT>;
      using ContributionT = typename AssemblerT::Contribution;

      System(const SystemDataT& system_data, const Model::StateData& state);
      ~System() override = default;

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

      int evaluateHessian(RealT        objective_factor,
                          const RealT* multipliers,
                          IdxT         multiplier_count) override;

      CsrMatrixT* getCsrHessian() const override;
      bool        hasHessian() override;

      const SystemDataT&      systemData() const;
      const Model::StateData& inputState() const;
      Model::StateData        solutionState() const;

    private:
      void validate() const;
      int  gatherVariables();
      void requireInitialized(const char* operation) const;

      SystemDataT      system_data_;
      Model::StateData input_state_;

      std::vector<std::unique_ptr<ComponentT>> components_;
      std::vector<ContributionT>               jacobian_contributions_;
      std::vector<ContributionT>               hessian_contributions_;
      std::vector<std::vector<RealT>>          local_multiplier_buffers_;

      VectorT variables_;
      VectorT variable_lower_bounds_;
      VectorT variable_upper_bounds_;
      VectorT objective_gradient_;
      VectorT constraints_;
      VectorT constraint_lower_bounds_;
      VectorT constraint_upper_bounds_;

      AssemblerT jacobian_assembler_;
      AssemblerT hessian_assembler_;

      IdxT  variable_count_{0};
      IdxT  constraint_count_{0};
      RealT objective_{0};
      bool  allocated_{false};
      bool  initialized_{false};
    };

    extern template class System<double, int>;
    extern template class System<double, long int>;
    extern template class System<double, std::size_t>;

  } // namespace OPF
} // namespace GridKit
