#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <sstream>
#include <string>

#include <GridKit/Model/OPF/System.hpp>
#include <GridKit/Model/OPF/SystemData.hpp>
#include <GridKit/Model/StateData.hpp>
#include <GridKit/Solver/Optimization/IpoptSolver.hpp>

namespace
{
  bool near(double value, double expected, double tolerance)
  {
    return std::abs(value - expected) <= tolerance;
  }

  GridKit::OPF::SystemData<> makeSystemData()
  {
    std::istringstream stream(R"(
      {
        "header": {
          "format_version": 0,
          "format_revision": 1,
          "case_name": "Two-bus Ipopt integration"
        },
        "params": {"freq_base": 60.0, "va_base": 100000000.0},
        "buses": [
          {
            "class": "Slack",
            "id": "B0",
            "params": {"number": 0, "kv": 230.0, "vmin": 1.0, "vmax": 1.0}
          },
          {
            "class": "Bus",
            "id": "B1",
            "params": {"number": 1, "kv": 230.0, "vmin": 0.95, "vmax": 1.05}
          }
        ],
        "devices": [
          {
            "class": "Branch",
            "id": "BR",
            "buses": {"from": 0, "to": 1},
            "params": {"R": 0.0, "X": 0.1, "G": 0.0, "B": 0.0, "smax": 2.0}
          },
          {
            "class": "Generator",
            "id": "G1",
            "buses": {"bus": 0},
            "params": {
              "pmin": 0.0, "pmax": 2.0,
              "qmin": 0.0, "qmax": 0.0,
              "mva": 100.0, "c2": 10.0
            }
          },
          {
            "class": "Generator",
            "id": "G2",
            "buses": {"bus": 0},
            "params": {
              "pmin": 0.0, "pmax": 2.0,
              "qmin": -2.0, "qmax": 2.0,
              "mva": 100.0, "c2": 20.0
            }
          },
          {
            "class": "Load",
            "id": "L",
            "buses": {"bus": 1},
            "params": {
              "pmin": -1.0, "pmax": -1.0,
              "qmin": 0.0, "qmax": 0.0
            }
          },
          {
            "class": "Shunt",
            "id": "SH",
            "buses": {"bus": 1},
            "params": {"G": 0.0, "B": 0.1}
          }
        ]
      })");
    return GridKit::OPF::parseSystemData(stream);
  }

  GridKit::Model::StateData makeStateData()
  {
    std::istringstream stream(R"(
      {
        "header": {"version": 1, "description": "Feasible OPF start"},
        "buses": {
          "bus_id_0": {
            "vr": 1.0,
            "vi": 0.0,
            "injections": {"preserved": {"ir": 3.0, "ii": 4.0}}
          },
          "bus_id_1": {
            "vr": 1.0001020304050608,
            "vi": -0.1,
            "injections": {}
          }
        },
        "devices": {
          "BR": {"open": false, "tap": 1.0, "phase": 0.0},
          "G1": {"online": true, "p": 0.5, "q": 0.0},
          "G2": {"online": true, "p": 0.5, "q": -0.0010203030301017},
          "L": {"online": true, "p": -1.0, "q": 0.0},
          "SH": {"online": true},
          "unrelated": {"active": false, "p": 7.0}
        }
      })");
    return GridKit::Model::parseStateData(stream);
  }
} // namespace

int main()
{
  using SystemT = GridKit::OPF::System<double, std::size_t>;

  SystemT system(makeSystemData(), makeStateData());
  if (system.allocate() != 0 || system.initialize() != 0)
  {
    std::cerr << "Could not allocate and initialize the OPF integration model\n";
    return 1;
  }

  if (system.size() != 8 || system.sizeResidual() != 6 || system.nnz() != 28)
  {
    std::cerr << "Unexpected OPF integration-model dimensions\n";
    return 1;
  }

  auto* jacobian = system.getCsrJacobian();
  if (jacobian == nullptr || jacobian->getNumRows() != 6
      || jacobian->getNumColumns() != 8 || jacobian->getNnz() != 28)
  {
    std::cerr << "Unexpected OPF integration-model Jacobian shape\n";
    return 1;
  }

  const auto* rows = jacobian->getRowData();
  const auto* cols = jacobian->getColData();
  for (std::size_t row = 0; row < system.sizeResidual(); ++row)
  {
    if (rows[row] > rows[row + 1])
    {
      std::cerr << "OPF Jacobian row pointers are not monotone\n";
      return 1;
    }
    for (std::size_t entry = rows[row] + 1; entry < rows[row + 1]; ++entry)
    {
      if (cols[entry - 1] >= cols[entry])
      {
        std::cerr << "OPF Jacobian columns are not sorted and unique\n";
        return 1;
      }
    }
  }

  AnalysisManager::IpoptInterface::Options<double, std::size_t> options;
  options.tolerance          = 1.0e-9;
  options.maximum_iterations = 200;
  options.print_level        = 0;

  AnalysisManager::IpoptInterface::Solver<double, std::size_t> solver(&system);
  const auto                                                   result = solver.solve(options);
  if (!result.solved())
  {
    std::cerr << "Ipopt did not solve the OPF integration model\n";
    return 1;
  }
  if (!std::isfinite(result.objective)
      || !std::isfinite(result.constraint_violation)
      || result.constraint_violation > 1.0e-7)
  {
    std::cerr << "Ipopt returned excessive OPF constraint violation\n";
    return 1;
  }

  if (system.evaluateResidual() != 0 || system.evaluateObjective() != 0)
  {
    std::cerr << "Could not reevaluate the accepted OPF point\n";
    return 1;
  }

  const auto* residual = system.getResidual().getData();
  double      max_residual{};
  bool        residuals_finite = true;
  for (std::size_t row = 0; row < system.sizeResidual(); ++row)
  {
    residuals_finite = residuals_finite && std::isfinite(residual[row]);
  }
  for (std::size_t row = 0; row < 4; ++row)
  {
    max_residual = std::max(max_residual, std::abs(residual[row]));
  }
  if (!residuals_finite || !std::isfinite(system.objective())
      || max_residual > 1.0e-7 || residual[4] < 0.0 || residual[4] > 4.0
      || residual[5] < 0.0 || residual[5] > 4.0
      || !near(system.objective(), 20.0 / 3.0, 1.0e-6))
  {
    std::cerr << "Accepted OPF point does not satisfy the analytic solution: objective="
              << system.objective() << ", max residual=" << max_residual << '\n';
    return 1;
  }

  const auto solution = system.solutionState();
  if (!near(*solution.devices.at("G1").p, 2.0 / 3.0, 1.0e-6)
      || !near(*solution.devices.at("G1").q, 0.0, 1.0e-6)
      || !near(*solution.devices.at("G2").p, 1.0 / 3.0, 1.0e-6)
      || !near(*solution.devices.at("G2").q, -0.00102030303030303, 1.0e-6)
      || !near(*solution.devices.at("L").p, -1.0, 1.0e-12)
      || solution.devices.at("BR").open != false
      || solution.devices.at("BR").tap != 1.0
      || solution.devices.at("BR").phase != 0.0
      || solution.devices.at("unrelated").active != false
      || solution.buses.at("bus_id_0").injections.count("G1") != 1
      || solution.buses.at("bus_id_0").injections.count("G2") != 1
      || solution.buses.at("bus_id_1").injections.count("L") != 1
      || solution.buses.at("bus_id_1").injections.count("SH") != 1
      || solution.buses.at("bus_id_0").injections.count("preserved") != 1)
  {
    std::cerr << "OPF solution-state output is incomplete or destructive\n";
    return 1;
  }

  const double vm1 = std::hypot(*solution.buses.at("bus_id_1").vr,
                                *solution.buses.at("bus_id_1").vi);
  const double va1 = std::atan2(*solution.buses.at("bus_id_1").vi,
                                *solution.buses.at("bus_id_1").vr);
  if (!near(vm1, 1.0050890861, 1.0e-6)
      || !near(va1, -0.0996585515, 1.0e-6))
  {
    std::cerr << "OPF solution voltage does not match the analytic high-voltage point\n";
    return 1;
  }

  std::ostringstream serialized;
  GridKit::Model::writeStateData(solution, serialized);
  std::istringstream serialized_input(serialized.str());
  const auto         round_trip = GridKit::Model::parseStateData(serialized_input);
  if (round_trip.devices.size() != solution.devices.size()
      || round_trip.buses.size() != solution.buses.size()
      || round_trip.devices.at("G1").p != solution.devices.at("G1").p
      || round_trip.buses.at("bus_id_1").injections.at("SH").ii
             != solution.buses.at("bus_id_1").injections.at("SH").ii)
  {
    std::cerr << "OPF solution state did not survive JSON round trip\n";
    return 1;
  }

  return 0;
}
