/**
 * @file EconomicDispatch.cpp
 * @brief AC optimal power flow of a MATPOWER case.
 */

#include <cmath>
#include <exception>
#include <format>
#include <string>

#include <IpIpoptApplication.hpp>

#include <GridKit/Model/OptimalPowerFlow/SystemModel.hpp>
#include <GridKit/Solver/Optimization/SteadyState.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace OPF = GridKit::OptimalPowerFlow;

using Log = GridKit::Utilities::Logger;

/// Relative tolerance of the optimal cost against a reference cost
constexpr double COST_TOLERANCE = 1.0e-8;

int runApplication(int argc, const char* argv[])
{
  Log::raiseVerbosity(Log::SUMMARY);

  if (argc < 2)
  {
    Log::error() << "No input file provided\n"
                 << "\n"
                 << "Usage:\n"
                 << "       EconomicDispatch <matpower-case-file> [reference-cost]\n"
                 << "\n";
    return 1;
  }

  OPF::SystemModel<double, size_t> model(OPF::parseMatpowerCase(argv[1]));
  if (model.allocate() != 0)
  {
    return 1;
  }

  // Exact sparse Jacobians and Hessians come from the model
  Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
  if (app->Initialize() != Ipopt::Solve_Succeeded)
  {
    return 1;
  }

  Ipopt::SmartPtr<Ipopt::TNLP> problem = new AnalysisManager::IpoptInterface::SteadyState<double, size_t>(&model);
  if (app->OptimizeTNLP(problem) != Ipopt::Solve_Succeeded)
  {
    Log::error() << "Ipopt did not find an optimal solution\n";
    return 1;
  }

  const double cost = model.objective();
  Log::summary() << std::format("Optimal cost {:.2f} $/h\n", cost);

  // A reference cost turns the study into a check
  if (argc > 2)
  {
    const double reference = std::stod(argv[2]);
    const double error     = std::abs(cost - reference) / std::abs(reference);
    if (error > COST_TOLERANCE)
    {
      Log::error() << std::format("Optimal cost differs from the reference {} $/h by {:.1e}\n", reference, error);
      return 1;
    }
  }

  return 0;
}

int main(int argc, const char* argv[])
{
  try
  {
    return runApplication(argc, argv);
  }
  catch (const std::exception& error)
  {
    Log::error() << "EconomicDispatch failed: " << error.what() << '\n';
  }

  return 1;
}
