#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include <IpIpoptApplication.hpp>

#include <GridKit/Model/OPF/System.hpp>
#include <GridKit/Model/OPF/SystemData.hpp>
#include <GridKit/Model/StateData.hpp>
#include <GridKit/Solver/Optimization/Ipopt/IpoptAdapter.hpp>

namespace
{
  bool near(double value, double expected, double tolerance)
  {
    return std::abs(value - expected) <= tolerance;
  }

  int fail(const std::string& message)
  {
    std::cerr << message << '\n';
    return 1;
  }

  GridKit::OPF::SystemData<> makeSystemData()
  {
    std::istringstream stream(R"(
      {
        "header": {
          "format_version": 0,
          "format_revision": 1,
          "case_name": "Two-bus exact-Hessian integration"
        },
        "params": {"freq_base": 60.0, "va_base": 100000000.0},
        "buses": [
          {
            "class": "Bus",
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
              "pmin": 0.0,
              "pmax": 2.0,
              "qmin": 0.0,
              "qmax": 0.0,
              "mva": 100.0,
              "c2": 10.0
            }
          },
          {
            "class": "Generator",
            "id": "G2",
            "buses": {"bus": 0},
            "params": {
              "pmin": 0.0,
              "pmax": 2.0,
              "qmin": -2.0,
              "qmax": 2.0,
              "mva": 100.0,
              "c2": 20.0
            }
          },
          {
            "class": "Load",
            "id": "L",
            "buses": {"bus": 1},
            "params": {
              "pmin": -1.0,
              "pmax": -1.0,
              "qmin": 0.0,
              "qmax": 0.0
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

  bool sameInjection(const GridKit::Model::InjectionState& left,
                     const GridKit::Model::InjectionState& right)
  {
    return left.ir == right.ir && left.ii == right.ii;
  }

  bool sameBus(const GridKit::Model::BusState& left,
               const GridKit::Model::BusState& right)
  {
    if (left.vr != right.vr || left.vi != right.vi
        || left.injections.size() != right.injections.size())
    {
      return false;
    }
    for (const auto& [id, injection] : left.injections)
    {
      const auto entry = right.injections.find(id);
      if (entry == right.injections.end()
          || !sameInjection(injection, entry->second))
      {
        return false;
      }
    }
    return true;
  }

  bool sameDevice(const GridKit::Model::DeviceState& left,
                  const GridKit::Model::DeviceState& right)
  {
    return left.active == right.active
           && left.online == right.online
           && left.open == right.open
           && left.p == right.p
           && left.q == right.q
           && left.tap == right.tap
           && left.phase == right.phase;
  }

  bool sameState(const GridKit::Model::StateData& left,
                 const GridKit::Model::StateData& right)
  {
    if (left.header.has_value() != right.header.has_value()
        || left.buses.size() != right.buses.size()
        || left.devices.size() != right.devices.size())
    {
      return false;
    }
    if (left.header
        && (left.header->version != right.header->version
            || left.header->time != right.header->time
            || left.header->created != right.header->created
            || left.header->description != right.header->description))
    {
      return false;
    }
    for (const auto& [id, bus] : left.buses)
    {
      const auto entry = right.buses.find(id);
      if (entry == right.buses.end() || !sameBus(bus, entry->second))
      {
        return false;
      }
    }
    for (const auto& [id, device] : left.devices)
    {
      const auto entry = right.devices.find(id);
      if (entry == right.devices.end()
          || !sameDevice(device, entry->second))
      {
        return false;
      }
    }
    return true;
  }

  template <typename VectorT>
  bool finiteVector(const VectorT& vector)
  {
    const auto* values = vector.getData();
    if (values == nullptr && vector.getSize() != 0)
    {
      return false;
    }
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
    if (values == nullptr && matrix->getNnz() != 0)
    {
      return false;
    }
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
} // namespace

int main()
{
  using SystemT  = GridKit::OPF::System<double, std::size_t>;
  using AdapterT = GridKit::Optimization::IpoptAdapter<double, std::size_t>;

  SystemT system(makeSystemData(), makeStateData());
  if (system.allocate() != 0 || system.initialize() != 0)
  {
    return fail("Could not allocate and initialize the OPF integration model");
  }
  if (system.size() != 8 || system.sizeConstraints() != 6
      || !system.hasJacobian() || !system.hasHessian()
      || system.getCsrJacobian() == nullptr
      || system.getCsrHessian() == nullptr)
  {
    return fail("Unexpected OPF integration-model derivative structure");
  }

  Ipopt::SmartPtr<Ipopt::IpoptApplication> application =
      IpoptApplicationFactory();
  if (!application->Options()->SetStringValue("derivative_test",
                                              "second-order")
      || !application->Options()->SetStringValue("derivative_test_print_all",
                                                 "yes")
      || !application->Options()->SetNumericValue(
          "derivative_test_perturbation", 1.0e-8)
      || !application->Options()->SetNumericValue("derivative_test_tol",
                                                  1.0e-6)
      || !application->Options()->SetNumericValue("tol", 1.0e-10)
      || !application->Options()->SetIntegerValue("max_iter", 200)
      || !application->Options()->SetIntegerValue("print_level", 5)
      || !application->Options()->SetIntegerValue("mumps_print_level", 0))
  {
    return fail("Could not configure exact-Hessian Ipopt integration options");
  }

  auto status = application->Initialize();
  if (status != Ipopt::Solve_Succeeded)
  {
    return fail("Could not initialize Ipopt for the OPF integration model");
  }

  Ipopt::SmartPtr<Ipopt::TNLP> problem = new AdapterT(system);
  status                               = application->OptimizeTNLP(problem);
  if (status != Ipopt::Solve_Succeeded)
  {
    return fail("Ipopt did not solve the exact-Hessian OPF integration model");
  }

  std::vector<double> multipliers(
      static_cast<std::size_t>(system.sizeConstraints()));
  for (std::size_t entry = 0; entry < multipliers.size(); ++entry)
  {
    multipliers[entry] = 0.125 * static_cast<double>(entry + 1);
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
    return fail("Could not reevaluate finite exact derivatives at the OPF solution");
  }

  const auto* constraints = system.constraints().getData();
  double      max_balance = 0.0;
  for (std::size_t row = 0; row < 4; ++row)
  {
    max_balance = std::max(max_balance, std::abs(constraints[row]));
  }
  if (max_balance > 1.0e-7
      || constraints[4] < -1.0e-10 || constraints[4] > 4.0 + 1.0e-10
      || constraints[5] < -1.0e-10 || constraints[5] > 4.0 + 1.0e-10
      || !near(system.objective(), 20.0 / 3.0, 1.0e-6))
  {
    return fail("The accepted OPF point does not satisfy the analytic solution");
  }

  const auto solution = system.solutionState();
  if (!solution.devices.at("G1").p
      || !solution.devices.at("G1").q
      || !solution.devices.at("G2").p
      || !solution.devices.at("G2").q
      || !near(*solution.devices.at("G1").p, 2.0 / 3.0, 1.0e-6)
      || !near(*solution.devices.at("G1").q, 0.0, 1.0e-6)
      || !near(*solution.devices.at("G2").p, 1.0 / 3.0, 1.0e-6)
      || !near(*solution.devices.at("G2").q,
               -0.00102030303030303,
               1.0e-6)
      || solution.devices.at("L").p != std::optional<double>{-1.0}
      || solution.devices.at("BR").open != std::optional<bool>{false}
      || solution.devices.at("BR").tap != std::optional<double>{1.0}
      || solution.devices.at("BR").phase != std::optional<double>{0.0}
      || solution.devices.at("unrelated").active
             != std::optional<bool>{false}
      || solution.devices.at("unrelated").p != std::optional<double>{7.0}
      || solution.buses.at("bus_id_0").injections.count("G1") != 1
      || solution.buses.at("bus_id_0").injections.count("G2") != 1
      || solution.buses.at("bus_id_1").injections.count("L") != 1
      || solution.buses.at("bus_id_1").injections.count("SH") != 1
      || solution.buses.at("bus_id_0").injections.count("preserved") != 1
      || solution.buses.at("bus_id_0").injections.at("preserved").ir
             != std::optional<double>{3.0}
      || solution.buses.at("bus_id_0").injections.at("preserved").ii
             != std::optional<double>{4.0})
  {
    return fail("OPF solution-state output is incomplete or destructive");
  }

  const auto& bus1 = solution.buses.at("bus_id_1");
  if (!bus1.vr || !bus1.vi)
  {
    return fail("OPF solution-state output is missing the solved bus voltage");
  }
  const double voltage_magnitude = std::hypot(*bus1.vr, *bus1.vi);
  const double voltage_angle     = std::atan2(*bus1.vi, *bus1.vr);
  if (!near(voltage_magnitude, 1.0050890861, 1.0e-6)
      || !near(voltage_angle, -0.0996585515, 1.0e-6))
  {
    return fail(
        "OPF solution voltage does not match the analytic high-voltage point");
  }

  std::ostringstream serialized;
  GridKit::Model::writeStateData(solution, serialized);
  std::istringstream serialized_input(serialized.str());
  const auto         round_trip = GridKit::Model::parseStateData(serialized_input);
  if (!sameState(solution, round_trip))
  {
    return fail(
        "The complete recognized OPF solution state did not survive JSON round trip");
  }

  return 0;
}
