/**
 * @file solve_pf.cpp
 *
 * Reads a MATPOWER v2 .m file passed as argv[1], runs the GridKit
 * Newton/KINSOL power flow solver, and prints per-bus results to stdout
 * in a simple machine-parseable format:
 *
 *   bus <bus_i>  V=<pu>  theta_deg=<deg>  type=<1|2|3>
 *
 * Exit code 0 = converged, non-zero = solver failure.
 *
 * Usage:
 *   solve_pf <path/to/case.m> [--output-m <path/to/solved.m>]
 *
 * If --output-m is given and the solve converges, a new .m file is written
 * with Vm/Va in mpc.bus replaced by the PF solution. Everything else
 * (mpc.gen, mpc.branch, gencost, comments) is copied through unchanged.
 */

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>
#include <iostream>
#include <map>
#include <string>
#include <utility>

#include <GridKit/Model/PowerFlow/Branch/Branch.hpp>
#include <GridKit/Model/PowerFlow/Bus/BusFactory.hpp>
#include <GridKit/Model/PowerFlow/Generator/GeneratorFactory.hpp>
#include <GridKit/Model/PowerFlow/Load/Load.hpp>
#include <GridKit/Model/PowerFlow/PowerFlowData.hpp>
#include <GridKit/Model/PowerFlow/SystemModelPowerFlow.hpp>
#include <GridKit/Solver/SteadyState/Kinsol.hpp>
#include "matpower_parser.hpp"  // local parser: handles underscore fields + extra columns

using namespace GridKit;
using namespace GridKit::PowerFlowData;
using namespace AnalysisManager::Sundials;
using namespace AnalysisManager;

int main(int argc, char* argv[])
{
  if (argc < 2) {
    std::cerr << "Usage: solve_pf <path/to/case.m> [--output-m <path>] [--flat-start]\n";
    return 1;
  }

  std::string m_path = argv[1];
  std::string out_m_path;
  bool flat_start = false;
  for (int i = 2; i < argc; ++i) {
    if (std::string(argv[i]) == "--output-m" && i + 1 < argc) {
      out_m_path = argv[++i];
    } else if (std::string(argv[i]) == "--flat-start") {
      flat_start = true;
    }
  }
  std::cerr << "Reading: " << m_path << "\n";

  // --- parse ---
  // Use the local parser (matpower_parser.hpp) rather than GridKit's
  // MatpowerParser.hpp to avoid modifying GridKit source. The local parser
  // handles extended MATPOWER formats (underscore field names, extra columns).
  SystemModelData<double, size_t> mp;
  try {
    UqPfSolver::readMatPowerFile(mp, m_path);
  } catch (const std::exception& e) {
    std::cerr << "Parse error: " << e.what() << "\n";
    return 2;
  }
  std::cerr << "Parsed: " << mp.bus.size() << " buses, "
            << mp.gen.size() << " gens, "
            << mp.branch.size() << " branches\n";

  // Flat start: override parsed Vm/Va with Vm=1.0, Va=0.0 for all buses.
  if (flat_start) {
    for (auto& bd : mp.bus) {
      bd.Vm = 1.0;
      bd.Va = 0.0;  // already stored in radians by local parser
    }
    std::cerr << "Flat start: Vm=1.0 Va=0.0 for all buses\n";
  }

  // --- assemble ---
  auto* sys = new SystemSteadyStateModel<double, size_t>(mp);
  sys->allocate();
  sys->initialize();
  std::cerr << "Model size (DOF): " << sys->size() << "\n";

  // --- solve ---
  auto* kinsol = new Kinsol<double, size_t>(sys);
  kinsol->configureSimulation();
  kinsol->getDefaultInitialCondition();   // warm start from parsed Vm/Va (deg->rad)
  int ret = 0;
  try {
    ret = kinsol->runSimulation();
  } catch (const std::exception& e) {
    std::cerr << "KINSOL exception: " << e.what() << "\n";
    ret = -99;
  }

  std::cerr << "KINSOL return code: " << ret << "\n";

  // Evaluate residual norm to decide convergence independently of KINSOL's
  // return code. KINSOL -5 (KIN_LINESEARCH_NONCONV) fires when warm-started
  // at the converged point -- residual is tiny but KINSOL throws. We accept
  // the result if ||f|| < tol regardless of the return code.
  sys->evaluateResidual();
  const auto& fvec = sys->getResidual();
  double fnorm = 0.0;
  for (size_t i = 0; i < fvec.size(); ++i)
    fnorm += fvec[i] * fvec[i];
  fnorm = std::sqrt(fnorm);
  std::cerr << "Residual 2-norm: " << fnorm << "\n";

  const double CONV_TOL = 1e-4;
  bool converged = (fnorm < CONV_TOL);
  std::cerr << (converged ? "CONVERGED" : "NOT CONVERGED")
            << "  (tol=" << CONV_TOL << ")\n";

  kinsol->printFinalStats();

  // --- print per-bus results (machine-parseable on stdout) ---
  // Only print if converged; otherwise the state is meaningless.
  if (converged) {
    std::map<size_t, std::pair<double, double>> bus_solution;
    for (const auto& bd : mp.bus) {
      auto* bus = sys->getBus(bd.bus_i);
      double v   = bus->V();
      double deg = bus->theta() * 180.0 / M_PI;
      bus_solution[bd.bus_i] = {v, deg};
      printf("bus %zu  V=%.8f  theta_deg=%.8f  type=%d\n",
             (size_t)bd.bus_i, v, deg, (int)bd.type);
    }

    // --- write solved .m if --output-m was given ---
    if (!out_m_path.empty()) {
      try {
        UqPfSolver::write_solved_m(m_path, out_m_path, bus_solution);
      } catch (const std::exception& e) {
        std::cerr << "write_solved_m error: " << e.what() << "\n";
      }
    }
  }

  delete kinsol;
  delete sys;
  return converged ? 0 : 1;
}
