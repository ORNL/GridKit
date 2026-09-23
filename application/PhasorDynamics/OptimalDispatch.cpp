/**
 * @file OptimalDispatch.cpp
 * @brief Optimal dispatch of a PhasorDynamics case, written as its initial
 * state.
 */

#include <exception>
#include <iostream>
#include <ranges>
#include <stdexcept>
#include <string>

#include <IpIpoptApplication.hpp>
#include <magic_enum/magic_enum.hpp>

#include <GridKit/Model/OptimalPowerFlow/SystemModel.hpp>
#include <GridKit/Model/OptimalPowerFlow/SystemModelData.hpp>
#include <GridKit/Model/StateData.hpp>
#include <GridKit/Solver/Optimization/OptimizationProblem.hpp>

#include "AnalysisUtilities.hpp"

namespace OPF = GridKit::OptimalPowerFlow;
namespace PD  = GridKit::PhasorDynamics;

using PD::json;
using PD::Log;

/**
 * @brief Dispatchable generators from PhasorDynamics machines
 */
template <typename MachineDataT>
void addGenerators(OPF::SystemModelData<>& data, const std::vector<MachineDataT>& machines)
{
  using Buses = typename MachineDataT::Buses;

  for (const auto& machine : machines)
  {
    auto& generator                           = data.generator.emplace_back();
    generator.id                              = machine.disambiguation_string;
    generator.buses[OPF::GeneratorBuses::bus] = machine.buses.at(Buses::bus);
  }
}

/**
 * @brief Optimal power flow network of a PhasorDynamics case
 *
 * Machines become generators, `LoadZIP` devices fixed loads, and `LoadZ`
 * devices shunts. Branches keep the PhasorDynamics parameters.
 */
OPF::SystemModelData<> network(const PD::SystemModelData<>& grid)
{
  using BusType = PD::SystemModelData<>::BusDataT::BusType;

  OPF::SystemModelData<> data;
  data.va_base = grid.va_base;

  for (const auto& bus : grid.bus)
  {
    auto& record    = data.bus.emplace_back();
    record.number   = bus.bus_id;
    record.infinite = bus.bus_type == BusType::SLACK;
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
      record.parameters[parameter.value()] = PD::realParameter(branch, key);
    }
  }

  addGenerators(data, grid.genrou);
  addGenerators(data, grid.gensal);
  addGenerators(data, grid.genclassical);
  addGenerators(data, grid.regca);

  for (const auto& load : grid.loadzip)
  {
    auto& record                      = data.load.emplace_back();
    record.id                         = load.disambiguation_string;
    record.buses[OPF::LoadBuses::bus] = load.buses.at(PD::LoadZIPBuses::bus);
  }

  // Y = 1 / (R + jX)
  for (const auto& load : grid.loadz)
  {
    const double r = PD::realParameter(load, PD::LoadZParameters::R);
    const double x = PD::realParameter(load, PD::LoadZParameters::X);

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
  PD::checkCommandLine(argc, "OptimalDispatch");
  const PD::StudyData study = PD::parseStudyData(argv[1]);
  if (study.dispatch_file.empty() || study.output_state_file.empty())
  {
    Log::error() << "OptimalDispatch requires dispatch_file and output_state_file\n";
    return 1;
  }

  OPF::SystemModelData<> data = network(study.model_data);
  OPF::applyMatpowerData(data, OPF::parseMatpowerData(study.dispatch_file));

  OPF::SystemModel<double, size_t> model(data, PD::extractState(study.model_data));
  if (model.allocate() != 0)
  {
    return 1;
  }

  // PhasorDynamics initialization rejects limits that relaxed bounds violate
  Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
  app->Options()->SetNumericValue("bound_relax_factor", 0.0);
  setOptions(*app, study.ipopt);
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

  GridKit::Model::writeStateData(model.solutionState(), study.output_state_file);
  std::cout << "State written to " << study.output_state_file << "\n";

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
