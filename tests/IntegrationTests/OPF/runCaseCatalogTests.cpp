#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include <IpIpoptApplication.hpp>

#include <GridKit/Model/OPF/System.hpp>
#include <GridKit/Model/OPF/SystemData.hpp>
#include <GridKit/Model/StateData.hpp>
#include <GridKit/Solver/Optimization/Ipopt/IpoptAdapter.hpp>

namespace
{
  using SystemT  = GridKit::OPF::System<double, std::size_t>;
  using AdapterT = GridKit::Optimization::IpoptAdapter<double, std::size_t>;

  template <typename VectorT>
  bool finiteVector(const VectorT& vector)
  {
    const auto* values = vector.getData();
    for (std::size_t entry = 0;
         entry < static_cast<std::size_t>(vector.getSize());
         ++entry)
    {
      if (!std::isfinite(values[entry]))
      {
        return false;
      }
    }
    return true;
  }

  template <typename MatrixT>
  bool finiteMatrix(MatrixT* matrix)
  {
    if (matrix == nullptr)
    {
      return false;
    }
    const auto* values = matrix->getValues();
    for (std::size_t entry = 0;
         entry < static_cast<std::size_t>(matrix->getNnz());
         ++entry)
    {
      if (!std::isfinite(values[entry]))
      {
        return false;
      }
    }
    return true;
  }

  double constraintViolation(const SystemT& system)
  {
    const auto& constraints = system.constraints();
    const auto& lower       = system.constraintLowerBounds();
    const auto& upper       = system.constraintUpperBounds();
    double      violation{};

    for (std::size_t row = 0;
         row < static_cast<std::size_t>(constraints.getSize());
         ++row)
    {
      violation = std::max(violation,
                           lower.getData()[row] - constraints.getData()[row]);
      violation = std::max(violation,
                           constraints.getData()[row] - upper.getData()[row]);
    }
    return violation;
  }

  int evaluate(SystemT& system)
  {
    std::vector<double> multipliers(
        static_cast<std::size_t>(system.sizeConstraints()));
    for (std::size_t entry = 0; entry < multipliers.size(); ++entry)
    {
      multipliers[entry] = 0.01 * static_cast<double>(entry % 17 + 1);
    }

    if (system.evaluateObjective() != 0
        || system.evaluateObjectiveGradient() != 0
        || system.evaluateConstraints() != 0
        || system.evaluateJacobian() != 0
        || system.evaluateHessian(0.75,
                                  multipliers.data(),
                                  system.sizeConstraints())
               != 0
        || !std::isfinite(system.objective())
        || !finiteVector(system.variables())
        || !finiteVector(system.objectiveGradient())
        || !finiteVector(system.constraints())
        || !finiteMatrix(system.getCsrJacobian())
        || !finiteMatrix(system.getCsrHessian()))
    {
      std::cerr << "Could not evaluate finite exact OPF derivatives\n";
      return 1;
    }
    return 0;
  }

  int solve(SystemT& system)
  {
    Ipopt::SmartPtr<Ipopt::IpoptApplication> application =
        IpoptApplicationFactory();
    application->Options()->SetNumericValue("tol", 1.0e-8);
    application->Options()->SetNumericValue("acceptable_tol", 1.0e-6);
    application->Options()->SetIntegerValue("max_iter", 500);
    application->Options()->SetIntegerValue("print_level", 0);
    application->Options()->SetIntegerValue("mumps_print_level", 0);

    auto status = application->Initialize();
    if (status != Ipopt::Solve_Succeeded)
    {
      std::cerr << "Could not initialize Ipopt\n";
      return 1;
    }

    Ipopt::SmartPtr<Ipopt::TNLP> problem = new AdapterT(system);
    status                               = application->OptimizeTNLP(problem);
    if (status != Ipopt::Solve_Succeeded
        && status != Ipopt::Solved_To_Acceptable_Level)
    {
      std::cerr << "Ipopt failed with status " << status << '\n';
      return 1;
    }

    if (evaluate(system) != 0)
    {
      std::cerr << "Could not evaluate exact derivatives at the accepted solution\n";
      return 1;
    }

    const double violation = constraintViolation(system);
    if (!std::isfinite(violation) || violation > 1.0e-6)
    {
      std::cerr << "Accepted solution has constraint violation "
                << violation << '\n';
      return 1;
    }

    std::cout << " status=" << status
              << " objective=" << system.objective()
              << " violation=" << violation;
    return 0;
  }
} // namespace

int main(int argc, char** argv)
{
  if (argc < 3 || argc > 4)
  {
    std::cerr << "Usage: " << argv[0]
              << " CASE.opf.json CASE.state.json [solve]\n";
    return 1;
  }

  try
  {
    SystemT system(GridKit::OPF::parseSystemData(std::filesystem::path(argv[1])),
                   GridKit::Model::parseStateData(
                       std::filesystem::path(argv[2])));
    if (system.allocate() != 0 || system.initialize() != 0)
    {
      std::cerr << "Could not allocate and initialize the OPF case\n";
      return 1;
    }
    if (!system.hasJacobian() || !system.hasHessian()
        || evaluate(system) != 0)
    {
      return 1;
    }

    std::cout << std::filesystem::path(argv[1]).stem().string()
              << " n=" << system.size()
              << " m=" << system.sizeConstraints()
              << " jacobian_nnz=" << system.getCsrJacobian()->getNnz()
              << " hessian_nnz=" << system.getCsrHessian()->getNnz();

    if (argc == 4)
    {
      if (std::string(argv[3]) != "solve" || solve(system) != 0)
      {
        return 1;
      }
    }
    std::cout << '\n';
  }
  catch (const std::exception& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }

  return 0;
}
