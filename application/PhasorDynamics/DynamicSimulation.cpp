#include <chrono>
#include <exception>
#include <filesystem>
#include <fstream>
#include <memory>
#include <vector>

#include <GridKit/Model/PhasorDynamics/BusFault/BusFault.hpp>
#include <GridKit/Model/PhasorDynamics/PartitionData.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/VariableMonitorController.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Solver/Dynamic/SplittingStep.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

#include "AnalysisUtilities.hpp"

using namespace AnalysisManager::Sundials;
using namespace GridKit::PhasorDynamics;
using namespace GridKit::Testing;

using scalar_type = double;
using real_type   = double;
using index_type  = size_t;

using SystemT = SystemModel<scalar_type, index_type>;
using IdaT    = Ida<scalar_type, index_type>;

/// Apply the study's IDA options and configure the solver.
void configure(IdaT& ida, const StudyData& study)
{
  ida.setTolerance(study.rel_tol, study.abs_tol);
  ida.setFixedStep(study.dt_fixed);
  ida.setMaxSteps(study.max_steps);
  ida.setMaxOrder(study.max_order);
  ida.setConsistentICType(study.consistent_ic_type);
  ida.configureSimulation();
}

/// Returns the simulation time in seconds, excluding setup.
real_type runMonolithic(const StudyData& study)
{
  SystemT system(study.model_data);
  system.allocate();

  IdaT ida(&system);
  configure(ida, study);

  const auto start = std::chrono::steady_clock::now();
  runStudy(study, ida, [&](std::size_t fault, bool on)
           { system.getBusFault(fault)->setStatus(on); });
  const auto stop = std::chrono::steady_clock::now();

  system.stopMonitor();
  return std::chrono::duration<real_type>(stop - start).count();
}

/// Returns the simulation time in seconds, excluding setup.
real_type runPartitioned(const StudyData& study)
{
  const auto parts = partitionSystemModelData(study.model_data, parsePartitionData(study.partition->file));

  std::vector<std::unique_ptr<SystemT>> systems;
  std::vector<std::unique_ptr<IdaT>>    solvers;
  for (const auto& data : parts.partitions)
  {
    auto& system = *systems.emplace_back(std::make_unique<SystemT>(data));
    system.allocate();
    configure(*solvers.emplace_back(std::make_unique<IdaT>(&system)), study);
  }

  SplittingStep<scalar_type, index_type> splitting; // destroyed before the partitions it integrates
  const auto                             couplings = connectPartitions(systems);
  for (std::size_t p = 0; p < systems.size(); ++p)
  {
    splitting.addPartition(*solvers[p], couplings[p], solvers[p]->context());
  }
  splitting.setNumThreads(study.partition->threads);
  splitting.setFixedStep(study.partition->dt);
  splitting.setCouplingTolerance(study.partition->tol);
  splitting.setTolerance(study.rel_tol, study.abs_tol);

  // One output in the intact case's column order, through the case's sinks.
  real_type                                            time = 0.0;
  GridKit::Model::VariableMonitorController<real_type> output(time);
  for (const auto& bus : study.model_data.bus)
  {
    output.addMonitor(systems[parts.bus_partition.at(bus.bus_id)]->getBus(bus.bus_id)->getMonitor());
  }
  for (const auto& [p, id] : parts.components)
  {
    output.addMonitor(systems[p]->getComponent(id)->getMonitor());
  }
  for (const auto& sink : study.model_data.monitor_sink)
  {
    output.addSink(sink);
  }
  output.start();
  splitting.setOutput([&](real_type t)
                      {
    time = t;
    output.print(); });
  splitting.configureSimulation();

  const auto start = std::chrono::steady_clock::now();
  runStudy(study, splitting, [&](std::size_t fault, bool on)
           {
    const auto& [p, id] = parts.faults.at(fault);
    systems[p]->getBusFault(id)->setStatus(on); });
  const auto stop = std::chrono::steady_clock::now();

  output.stop();
  Log::summary() << "Splitting steps: " << splitting.numSteps() << " (rejected " << splitting.numRejectedSteps() << ")\n";
  for (std::size_t p = 0; p < systems.size(); ++p)
  {
    Log::summary() << "Partition " << p + 1 << ": " << systems[p]->size() << " variables, "
                   << splitting.partitionTime(p) << " seconds\n";
  }
  return std::chrono::duration<real_type>(stop - start).count();
}

int runApplication(int argc, const char* argv[])
{
  // Print summaries, such as the run time, without lowering a higher verbosity
  Log::raiseVerbosity(Log::SUMMARY);

  // Study file
  checkCommandLine(argc, "DynamicSimulation");
  auto study = parseStudyData(argv[1]);

  real_type elapsed = 0.0;
  if (study.partition)
  {
    elapsed = runPartitioned(study);
  }
  else
  {
    elapsed = runMonolithic(study);
  }

  // Generate aggregate errors comparing variable output to reference solution
  TestStatus status = checkErrors(study);

  // Report run time
  Log::summary() << "Complete in " << elapsed << " seconds\n";

  return status.get();
}

int main(int argc, const char* argv[])
{
  try
  {
    return runApplication(argc, argv);
  }
  catch (const std::exception& error)
  {
    Log::error() << "DynamicSimulation failed: " << error.what() << '\n';
  }

  return 1;
}
