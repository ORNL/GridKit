/**
 * @file OptimalDispatch.cpp
 * @brief Optimal dispatch of a PhasorDynamics case, written as its initial
 * state.
 */

#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <ranges>
#include <stdexcept>
#include <string>
#include <variant>

#include <IpIpoptApplication.hpp>
#include <magic_enum/magic_enum.hpp>
#include <nlohmann/json.hpp>

#include <GridKit/Model/OptimalPowerFlow/SystemModel.hpp>
#include <GridKit/Model/OptimalPowerFlow/SystemModelData.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Model/StateData.hpp>
#include <GridKit/Solver/Optimization/OptimizationProblem.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace OPF = GridKit::OptimalPowerFlow;
namespace PD  = GridKit::PhasorDynamics;
namespace fs  = std::filesystem;

using json = nlohmann::json;
using Log  = GridKit::Utilities::Logger;
using GridKit::Model::StateData;

/// Numeric case parameter, or its model default if omitted
template <typename DataT>
typename DataT::RealT realParameter(const DataT& data, typename DataT::Parameters key, typename DataT::RealT fallback = 0.0)
{
  const auto entry = data.parameters.find(key);
  if (entry == data.parameters.end())
  {
    return fallback;
  }
  if (const auto* value = std::get_if<typename DataT::RealT>(&entry->second))
  {
    return *value;
  }
  if (const auto* value = std::get_if<typename DataT::IdxT>(&entry->second))
  {
    return static_cast<typename DataT::RealT>(*value);
  }
  throw std::invalid_argument(data.disambiguation_string + ": parameter " + std::string(magic_enum::enum_name(key)) + " must be numeric");
}

/// Preserve a supplied current pair, or initialize it from the case power
void initializeCurrent(StateData& state, const std::string& id, size_t bus, double p, double q)
{
  const auto* record = state.device(id);
  const bool  has_ir = record != nullptr && record->values.contains("ir");
  const bool  has_ii = record != nullptr && record->values.contains("ii");
  if (has_ir != has_ii)
  {
    throw std::invalid_argument(id + ": state must supply both ir and ii");
  }
  if (!has_ir)
  {
    GridKit::Model::setTerminalCurrent(state, id, bus, 0, 1, p, q);
  }
}

/**
 * @brief Dispatchable generators from PhasorDynamics machines
 */
template <typename MachineDataT>
void addGenerators(OPF::SystemModelData<>& data, StateData& state, const std::vector<MachineDataT>& machines)
{
  using Buses      = typename MachineDataT::Buses;
  using Parameters = typename MachineDataT::Parameters;

  for (const auto& machine : machines)
  {
    auto& generator                           = data.generator.emplace_back();
    generator.id                              = machine.disambiguation_string;
    generator.buses[OPF::GeneratorBuses::bus] = machine.buses.at(Buses::bus);

    const auto* record = state.device(generator.id);
    if (record != nullptr && !record->flag("online", true))
    {
      throw std::invalid_argument("Offline machine " + generator.id + " is not supported yet");
    }
    initializeCurrent(state, generator.id, machine.buses.at(Buses::bus), realParameter(machine, Parameters::p0), realParameter(machine, Parameters::q0));
  }
}

/**
 * @brief Optimal power flow network of a PhasorDynamics case
 *
 * Machines become generators, `LoadZIP` devices fixed loads, and `LoadZ`
 * devices shunts. Branches keep the PhasorDynamics parameters. Missing state
 * voltages and current pairs are initialized from the case in the same pass.
 */
OPF::SystemModelData<> network(const PD::SystemModelData<>& grid, StateData& state)
{
  using BusType = PD::SystemModelData<>::BusDataT::BusType;

  OPF::SystemModelData<> data;
  data.va_base = grid.va_base;

  for (const auto& bus : grid.bus)
  {
    auto& record    = data.bus.emplace_back();
    record.number   = bus.bus_id;
    record.infinite = bus.bus_type == BusType::SLACK;

    auto& voltage = state.buses[GridKit::Model::busKey(bus.bus_id)].values;
    voltage.try_emplace("vr", bus.Vr0);
    voltage.try_emplace("vi", bus.Vi0);
  }

  for (const auto& branch : grid.branch)
  {
    auto& record                         = data.branch.emplace_back();
    record.id                            = branch.disambiguation_string;
    record.buses[OPF::BranchBuses::bus1] = branch.buses.at(PD::BranchBuses::bus1);
    record.buses[OPF::BranchBuses::bus2] = branch.buses.at(PD::BranchBuses::bus2);
    for (const auto& key : std::views::keys(branch.parameters))
    {
      const auto parameter                 = magic_enum::enum_cast<OPF::BranchParameters>(magic_enum::enum_name(key));
      record.parameters[parameter.value()] = realParameter(branch, key);
    }
  }

  addGenerators(data, state, grid.genrou);
  addGenerators(data, state, grid.gensal);
  addGenerators(data, state, grid.genclassical);
  addGenerators(data, state, grid.regca);

  for (const auto& load : grid.loadzip)
  {
    auto& record                      = data.load.emplace_back();
    record.id                         = load.disambiguation_string;
    record.buses[OPF::LoadBuses::bus] = load.buses.at(PD::LoadZIPBuses::bus);
    initializeCurrent(state, record.id, load.buses.at(PD::LoadZIPBuses::bus), -realParameter(load, PD::LoadZIPParameters::Pnom), -realParameter(load, PD::LoadZIPParameters::Qnom));
  }

  // Y = 1 / (R + jX)
  for (const auto& load : grid.loadz)
  {
    const double r = realParameter(load, PD::LoadZParameters::R, 0.1);
    const double x = realParameter(load, PD::LoadZParameters::X, 0.01);

    auto& record                               = data.shunt.emplace_back();
    record.id                                  = load.disambiguation_string;
    record.buses[OPF::ShuntBuses::bus]         = load.buses.at(PD::LoadZBuses::bus);
    record.parameters[OPF::ShuntParameters::G] = r / (r * r + x * x);
    record.parameters[OPF::ShuntParameters::B] = -x / (r * r + x * x);
  }

  return data;
}

/**
 * @brief Set Ipopt options of string, integer, or real type
 */
void setOptions(Ipopt::IpoptApplication& app, const json& options)
{
  for (const auto& [name, value] : options.items())
  {
    bool set = false;
    if (value.is_string())
    {
      set = app.Options()->SetStringValue(name, value.get<std::string>());
    }
    else if (value.is_number_integer())
    {
      set = app.Options()->SetIntegerValue(name, value.get<int>());
    }
    else if (value.is_number())
    {
      set = app.Options()->SetNumericValue(name, value.get<double>());
    }

    if (!set)
    {
      throw std::invalid_argument("Invalid Ipopt option " + name);
    }
  }
}

int runApplication(int argc, const char* argv[])
{
  if (argc < 2)
  {
    std::cerr << "Usage: OptimalDispatch <json-input-file>\n";
    return 1;
  }

  const fs::path file = argv[1];
  std::ifstream  stream(file);
  if (!stream)
  {
    throw std::runtime_error("Could not open dispatch study " + file.string());
  }
  const json options = json::parse(stream);
  const auto input   = [&](const char* key)
  {
    const auto path = options.at(key).get<fs::path>();
    return path.is_absolute() ? path : file.parent_path() / path;
  };
  const auto output = options.at("output_state_file").get<fs::path>();
  const auto grid   = PD::parseSystemModelData(input("system_model_file"));

  StateData state;
  if (options.contains("state_file") && !options.at("state_file").get<fs::path>().empty())
  {
    state = GridKit::Model::parseStateData(input("state_file"));
  }
  auto data = network(grid, state);
  OPF::applyMatpowerData(data, OPF::parseMatpowerData(input("dispatch_file")));

  OPF::SystemModel<double, size_t> model(data, state);
  if (model.allocate() != 0)
  {
    return 1;
  }

  // PhasorDynamics initialization rejects limits that relaxed bounds violate
  Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
  app->Options()->SetNumericValue("bound_relax_factor", 0.0);
  setOptions(*app, options.value("ipopt", json::object()));
  if (app->Initialize() != Ipopt::Solve_Succeeded)
  {
    Log::error() << "Ipopt initialization failed\n";
    return 1;
  }

  Ipopt::SmartPtr<Ipopt::TNLP> problem = new AnalysisManager::IpoptInterface::OptimizationProblem<double, size_t>(&model);

  const Ipopt::ApplicationReturnStatus status = app->OptimizeTNLP(problem);
  if (status != Ipopt::Solve_Succeeded && status != Ipopt::Solved_To_Acceptable_Level)
  {
    Log::error() << "Ipopt did not converge, status " << magic_enum::enum_name(status) << "\n";
    return 1;
  }

  model.evaluateObjective();
  std::cout << "\nOptimal cost " << model.objective() << "\n";

  GridKit::Model::writeStateData(model.solutionState(), output);
  std::cout << "State written to " << output << "\n";

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
    Log::error() << "OptimalDispatch failed: " << error.what() << '\n';
  }

  return 1;
}
