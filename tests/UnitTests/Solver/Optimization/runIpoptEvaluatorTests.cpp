#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>

#include <GridKit/MemoryUtilities/MemoryUtils.hpp>
#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Solver/Optimization/IpoptSolver.hpp>

namespace
{
  class Quadratic final : public GridKit::Model::Evaluator<double, std::size_t>
  {
  public:
    explicit Quadratic(bool constrained)
      : constrained_(constrained)
    {
    }

    int allocate() override
    {
      int status  = variables_.resize(1);
      status     |= constraints_.resize(sizeResidual());
      status     |= objective_gradient_.resize(1);

      if (!constrained_)
      {
        return status;
      }

      jacobian_  = std::make_unique<CsrMatrixT>(1, 1, 1);
      status    |= jacobian_->allocateMatrixData(GridKit::memory::HOST);
      if (status != 0)
      {
        return status;
      }

      jacobian_->getRowData()[0] = 0;
      jacobian_->getRowData()[1] = 1;
      jacobian_->getColData()[0] = 0;
      jacobian_->getValues()[0]  = 1.0;
      return jacobian_->setUpdated(GridKit::memory::HOST);
    }

    int initialize() override
    {
      variables_.getData()[0] = 4.0;
      return variables_.setDataUpdated();
    }

    int evaluateResidual() override
    {
      if (!constrained_)
      {
        return 0;
      }
      constraints_.getData()[0] = variables_.getData()[0];
      return constraints_.setDataUpdated();
    }

    int evaluateJacobian() override
    {
      objective_gradient_.getData()[0] = 2.0 * (variables_.getData()[0] - 1.0);
      int status                       = 0;
      if (constrained_)
      {
        jacobian_->getValues()[0] = 1.0;
        status                    = jacobian_->setUpdated(GridKit::memory::HOST);
      }
      status |= objective_gradient_.setDataUpdated();
      return status;
    }

    int evaluateObjective() override
    {
      const double delta = variables_.getData()[0] - 1.0;
      objective_         = delta * delta;
      return 0;
    }

    std::size_t size() override
    {
      return 1;
    }

    std::size_t sizeResidual() override
    {
      return constrained_ ? 1 : 0;
    }

    std::size_t nnz() override
    {
      return constrained_ ? 1 : 0;
    }

    bool hasJacobian() override
    {
      return constrained_ && jacobian_ != nullptr;
    }

    CsrMatrixT* getCsrJacobian() const override
    {
      return jacobian_.get();
    }

    bool hasObjective() const override
    {
      return true;
    }

    double objective() const override
    {
      return objective_;
    }

    VectorT& y() override
    {
      return variables_;
    }

    const VectorT& y() const override
    {
      return variables_;
    }

    VectorT& getResidual() override
    {
      return constraints_;
    }

    const VectorT& getResidual() const override
    {
      return constraints_;
    }

    VectorT& getObjectiveGradient() override
    {
      return objective_gradient_;
    }

    const VectorT& getObjectiveGradient() const override
    {
      return objective_gradient_;
    }

  private:
    VectorT                     variables_;
    VectorT                     constraints_;
    VectorT                     objective_gradient_;
    std::unique_ptr<CsrMatrixT> jacobian_;
    double                      objective_{};
    bool                        constrained_{false};
  };
} // namespace

int main()
{
  Quadratic constrained_model(true);
  if (constrained_model.allocate() != 0 || constrained_model.initialize() != 0)
  {
    std::cerr << "Could not initialize Ipopt evaluator test model\n";
    return 1;
  }

  AnalysisManager::IpoptInterface::Solver<double, std::size_t> constrained_solver(&constrained_model);
  const auto                                                   constrained_result = constrained_solver.solve();

  const bool correct_constrained_solution  = std::abs(constrained_model.y().getData()[0]) < 1.0e-8;
  const bool correct_constrained_objective = std::abs(constrained_result.objective - 1.0) < 1.0e-8;
  const bool constrained_result_sizes      = constrained_result.primal.size() == 1 && constrained_result.constraints.size() == 1 && constrained_result.constraint_duals.size() == 1 && constrained_result.lower_bound_duals.size() == 1 && constrained_result.upper_bound_duals.size() == 1;

  if (!constrained_result.solved() || !correct_constrained_solution || !correct_constrained_objective || !constrained_result_sizes)
  {
    std::cerr << "Ipopt evaluator constrained solve did not return the expected solution\n";
    return 1;
  }

  Quadratic unconstrained_model(false);
  if (unconstrained_model.allocate() != 0 || unconstrained_model.initialize() != 0)
  {
    std::cerr << "Could not initialize unconstrained Ipopt evaluator test model\n";
    return 1;
  }

  AnalysisManager::IpoptInterface::Solver<double, std::size_t> unconstrained_solver(&unconstrained_model);
  const auto                                                   unconstrained_result = unconstrained_solver.solve();

  const bool correct_unconstrained_solution =
      std::abs(unconstrained_model.y().getData()[0] - 1.0) < 1.0e-8;
  if (!unconstrained_result.solved() || !correct_unconstrained_solution || !unconstrained_result.constraints.empty())
  {
    std::cerr << "Ipopt evaluator unconstrained solve did not return the expected solution\n";
    return 1;
  }
  return 0;
}
