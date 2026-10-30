#include <chrono>
#include <exception>
#include <filesystem>
#include <future>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <GridKit/Model/PhasorDynamics/Branch/Branch.hpp>
#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/Gas/GasBusBase.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

#include "AnalysisUtilities.hpp"

using Clock = std::chrono::high_resolution_clock;
using Dur   = std::chrono::duration<double>;

using Log = GridKit::Utilities::Logger;

using namespace AnalysisManager::Sundials;
using namespace GridKit::PhasorDynamics;
using namespace GridKit::Testing;

using scalar_type = double;
using real_type   = double;
using index_type  = size_t;

TestStatus runStudy(StudyData study_data)
{
  // Instantiate system
  SystemModel<scalar_type, index_type> sys(study_data.model_data);
  sys.allocate();

  // Set up simulation
  Ida<scalar_type, index_type> ida(&sys);
  ida.setTolerance(study_data.rel_tol, study_data.abs_tol);
  ida.setFixedStep(study_data.dt_fixed);
  ida.setMaxSteps(study_data.max_steps);
  ida.setMaxOrder(study_data.max_order);
  ida.setConsistentICType(study_data.consistent_ic_type);
  ida.configureSimulation();

  using EventType = SystemEvent::Type;

  // Initialize simulation for first run
  real_type dt_monitor = study_data.dt_monitor;
  real_type final_time = study_data.tmax;
  ida.initializeSimulation(0.0, false);

  for (const auto& event : study_data.events)
  {
    // Run to event time
    ida.runSimulation(event.time, dt_monitor);

    // Set up run for event (to start at event time)
    int ret = 0;

    switch (event.type)
    {
    case EventType::FAULT_ON:
      ret = sys.getBus(event.bus)->setFault(true, event.R, event.X);
      break;
    case EventType::FAULT_OFF:
      ret = sys.getBus(event.bus)->clearFault();
      break;
    case EventType::IN_SERVICE:
      sys.getBranch(event.device)->setInService(true);
      break;
    case EventType::OUT_OF_SERVICE:
      sys.getBranch(event.device)->setInService(false);
      break;
    case EventType::WITHDRAWAL:
      ret = sys.getGasBus(event.bus)->setWithdrawal(event.w);
      break;
    }

    if (ret != 0)
    {
      Log::error() << "Event failed for bus " << event.bus
                   << " at time " << event.time << '\n';
      return {false};
    }

    // Re-initialize simulation at event time
    ida.initializeSimulation(event.time, true);
  }

  // Run to final time
  ida.runSimulation(final_time, dt_monitor);

  // Stop the variable monitor
  sys.stopMonitor();

  return checkErrors(study_data, false);
}

TestStatus singleFaultStudy(std::size_t fault_id, StudyData study_data)
{
  // Retarget fault events while preserving branch-switching events
  for (auto& event : study_data.events)
  {
    if (event.type == SystemEvent::Type::FAULT_ON || event.type == SystemEvent::Type::FAULT_OFF)
    {
      event.bus = study_data.model_data.bus[fault_id].bus_id;
    }
  }

  // Make distinct output files
  for (auto& sink : study_data.model_data.monitor_sink)
  {
    auto path      = std::filesystem::path(sink.file_name);
    auto ext       = path.extension().string();
    auto name      = path.stem().string();
    sink.file_name = name + "_" + std::to_string(fault_id) + ext;
  }

  try
  {
    return runStudy(study_data);
  }
  catch (...)
  {
    Log::warning() << "exception caught at fault id: " << fault_id << std::endl;
    return {false};
  }
}

void runStudySerial(const StudyData& study_data, std::vector<TestStatus>& stat_vec)
{
  for (std::size_t i = 0; i < study_data.model_data.bus.size(); ++i)
  {
    auto stat   = singleFaultStudy(i, study_data);
    stat_vec[i] = stat;
  }
}

#ifdef GRIDKIT_ENABLE_THREADS
void runStudyAsync(const StudyData& study_data, std::vector<TestStatus>& stat_vec)
{
  auto n_faults = study_data.model_data.bus.size();

  std::vector<std::future<TestStatus>> futures;
  futures.reserve(n_faults);
  for (std::size_t i = 0; i < n_faults; ++i)
  {
    futures.emplace_back(
        std::async(std::launch::async, singleFaultStudy, i, study_data));
  }

  for (std::size_t i = 0; i < n_faults; ++i)
  {
    auto stat   = futures[i].get();
    stat_vec[i] = stat;
  }
}
#endif

#ifdef _OPENMP
void runStudyOpenMP(const StudyData& study_data, std::vector<TestStatus>& stat_vec)
{
  auto n_faults = study_data.model_data.bus.size();
#pragma omp parallel for
  for (std::size_t i = 0; i < n_faults; ++i)
  {
    auto stat   = singleFaultStudy(i, study_data);
    stat_vec[i] = stat;
  }
}
#endif

int runApplication(int argc, const char* argv[])
{
  // Print summaries, such as the run time, without lowering a higher verbosity
  Log::raiseVerbosity(Log::SUMMARY);

  // Study file
  checkCommandLine(argc, "ContingencyAnalysis");
  auto study_data = parseStudyData(argv[1]);

  const auto start = Clock::now();

  auto faults   = study_data.model_data.bus;
  auto stat_vec = std::vector<TestStatus>(faults.size(), true);

  // Use std::async if threads are available (so far, std::async has out-performed OpenMP)
  // Otherwise, use OpenMP if available
  // Fall back to serial if neither threads or OpenMP are available
#if defined(GRIDKIT_ENABLE_THREADS)
  runStudyAsync(study_data, stat_vec);
#elif defined(_OPENMP)
  runStudyOpenMP(study_data, stat_vec);
#else
  runStudySerial(study_data, stat_vec);
#endif

  const auto stop = Clock::now();
  const auto dur  = std::chrono::duration<double>(stop - start);
  Log::summary() << "Complete in " << dur << "\n";

  TestStatus status;
  for (std::size_t i = 0; i < stat_vec.size(); ++i)
  {
    status *= stat_vec[i];
    if (!stat_vec[i])
    {
      Log::error() << "Study failed for fault: "
                   << faults[i].name << '\n';
    }
  }

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
    Log::error() << "ContingencyAnalysis failed: " << error.what() << '\n';
  }

  return 1;
}
