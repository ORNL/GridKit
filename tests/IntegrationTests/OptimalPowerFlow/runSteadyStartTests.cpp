/**
 * @file runSteadyStartTests.cpp
 * @brief PhasorDynamics starts in steady state from an optimal power flow state.
 *
 * Usage: test_opf_steady_start <case> <state> <residual tolerance> <drift tolerance>
 */

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Model/StateData.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace
{
  using SystemModelT = GridKit::PhasorDynamics::SystemModel<double, size_t>;

  /// Simulated time over which the state must hold [s]
  constexpr double STEADY_TIME = 1.0;

  /// Solver tolerances of `DynamicSimulation` by default
  constexpr double REL_TOL = 1.0e-7;
  constexpr double ABS_TOL = 1.0e-9;

  double maxResidual(SystemModelT& system)
  {
    system.evaluateResidual();

    const auto&   f      = system.getResidual();
    const double* values = f.getData();

    double norm = 0.0;
    for (size_t i = 0; i < f.getSize(); ++i)
    {
      norm = std::max(norm, std::abs(values[i]));
    }
    return norm;
  }

  double maxDrift(const std::vector<double>& y0, const double* y)
  {
    double drift = 0.0;
    for (size_t i = 0; i < y0.size(); ++i)
    {
      drift = std::max(drift, std::abs(y[i] - y0[i]) / (1.0 + std::abs(y0[i])));
    }
    return drift;
  }
} // namespace

int main(int argc, const char* argv[])
{
  using namespace GridKit;

  if (argc != 5)
  {
    std::cerr << "Usage: " << argv[0] << " <case> <state> <residual tolerance> <drift tolerance>\n";
    return 1;
  }

  const double residual_tol = std::atof(argv[3]);
  const double drift_tol    = std::atof(argv[4]);

  Testing::TestStatus success = "steadyStart";

  auto data = PhasorDynamics::parseSystemModelData(std::filesystem::path(argv[1]));
  PhasorDynamics::applyState(data, Model::parseStateData(std::filesystem::path(argv[2])));

  SystemModelT system(data);
  system.allocate();

  const double residual = maxResidual(system);
  std::cout << "Residual at t = 0: " << residual << "\n";
  success *= residual < residual_tol;

  const double*             y = system.y().getData();
  const std::vector<double> y0(y, y + system.size());

  AnalysisManager::Sundials::Ida<double, size_t> ida(&system);
  ida.setTolerance(REL_TOL, ABS_TOL);
  ida.configureSimulation();
  ida.initializeSimulation(0.0);
  ida.runSimulation(STEADY_TIME);

  const double drift = maxDrift(y0, system.y().getData());
  std::cout << "Relative drift after " << STEADY_TIME << " s: " << drift << "\n";
  success *= drift < drift_tol;

  return success.report();
}
