/**
 * @file EconomicDispatch.cpp
 * @brief AC optimal power flow of a JSON case.
 */

#include <cmath>
#include <exception>
#include <format>

#include <IpIpoptApplication.hpp>

#include <GridKit/Model/OptimalPowerFlow/SystemModel.hpp>
#include <GridKit/Solver/Optimization/IpoptInterface.hpp>
#include <GridKit/Utilities/CliArgs/CliArgs.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace OPF = GridKit::OptimalPowerFlow;

using GridKit::Utilities::ArgType;
using GridKit::Utilities::CliArgs;
using Log = GridKit::Utilities::Logger;

int runApplication(int argc, const char* argv[])
{
  Log::raiseVerbosity(Log::SUMMARY);

  CliArgs args{{.name = {"--case-file", "-c"}, .help = "JSON case file", .required = true},
               {.name = {"--reference-cost", "-r"}, .help = "Optimal cost to check [$/h]", .type = ArgType::Real},
               {.name     = {"--solver-tolerance"},
                .help     = "Ipopt convergence tolerance",
                .type     = ArgType::Real,
                .defaults = 1.0e-8},
               {.name     = {"--tolerance", "-t"},
                .help     = "Relative tolerance of the cost check",
                .type     = ArgType::Real,
                .defaults = 1.0e-6}};
  args.parseArgs(argc, argv);

  OPF::SystemModel<double, size_t> model(OPF::parseSystemModelData(args["case-file"]()));
  if (model.allocate() != 0 || model.initialize() != 0)
  {
    return 1;
  }

  // Exact sparse Jacobians and Hessians come from the model
  Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
  if (!app->Options()->SetNumericValue("tol", args["solver-tolerance"].as<double>()))
  {
    return 1;
  }
  if (app->Initialize() != Ipopt::Solve_Succeeded)
  {
    return 1;
  }

  Ipopt::SmartPtr<Ipopt::TNLP> problem = new AnalysisManager::SteadyState::IpoptInterface<double, size_t>(&model);
  if (app->OptimizeTNLP(problem) != Ipopt::Solve_Succeeded)
  {
    Log::error() << "Ipopt did not find an optimal solution\n";
    return 1;
  }

  const double cost = model.f();
  Log::summary() << std::format("Optimal cost {:.2f} $/h\n", cost);

  // A reference cost turns the study into a check
  if (!args["reference-cost"].empty())
  {
    const double reference = args["reference-cost"].as<double>();
    const double error     = std::abs(cost - reference) / std::abs(reference);
    if (error > args["tolerance"].as<double>())
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
