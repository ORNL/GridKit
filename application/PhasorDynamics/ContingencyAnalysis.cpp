#include <chrono>
#include <cmath>
#include <filesystem>
#include <future>

#ifdef _OPENMP
#include <omp.h>
#endif

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
  ida.configureSimulation();

  // Initilize simultation for first run
  real_type dt_monitor = study_data.dt_monitor;
  real_type final_time = study_data.tmax;
  ida.initializeSimulation(0.0, false);

  for (std::size_t i = 0; i < study_data.events.size();)
  {
    const real_type event_time = study_data.events[i].time;

    // Run to event time
    ida.runSimulation(event_time, dt_monitor);

    // Apply every input change scheduled for this time before reinitializing.
    do
    {
      const auto& event = study_data.events[i];
      sys.setInput(event.device_id, event.input, event.value);
      ++i;
    } while (i < study_data.events.size()
             && study_data.events[i].time == event_time);

    // Re-initialize simulation at event time
    ida.initializeSimulation(event_time, true);
  }

  // Run to final time
  ida.runSimulation(final_time, dt_monitor);

  // Stop the variable monitor
  sys.stopMonitor();

  return checkErrors(study_data, false);
}

TestStatus singleFaultStudy(const std::string& fault_id,
                            std::size_t        fault_index,
                            StudyData          study_data)
{
  // Retarget fault-activation events while preserving all other input events.
  for (auto& event : study_data.events)
  {
    if (event.input == "active")
    {
      event.device_id = fault_id;
    }
  }

  // Make distinct output files
  for (auto& sink : study_data.model_data.monitor_sink)
  {
    auto path      = std::filesystem::path(sink.file_name);
    auto ext       = path.extension().string();
    auto name      = path.stem().string();
    sink.file_name = name + "_" + std::to_string(fault_index) + ext;
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
  for (std::size_t i = 0; i < study_data.model_data.bus_fault.size(); ++i)
  {
    auto stat = singleFaultStudy(
        study_data.model_data.bus_fault[i].disambiguation_string,
        i,
        study_data);
    stat_vec[i] = stat;
  }
}

#ifdef GRIDKIT_ENABLE_THREADS
void runStudyAsync(const StudyData& study_data, std::vector<TestStatus>& stat_vec)
{
  auto n_faults = study_data.model_data.bus_fault.size();

  std::vector<std::future<TestStatus>> futures;
  futures.reserve(n_faults);
  for (std::size_t i = 0; i < n_faults; ++i)
  {
    futures.emplace_back(
        std::async(
            std::launch::async,
            singleFaultStudy,
            study_data.model_data.bus_fault[i].disambiguation_string,
            i,
            study_data));
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
  auto n_faults = study_data.model_data.bus_fault.size();
#pragma omp parallel for
  for (std::size_t i = 0; i < n_faults; ++i)
  {
    auto stat = singleFaultStudy(
        study_data.model_data.bus_fault[i].disambiguation_string,
        i,
        study_data);
    stat_vec[i] = stat;
  }
}
#endif

int main(int argc, const char* argv[])
{
  // Study file
  checkCommandLine(argc, "ContingencyAnalysis");
  auto study_data = parseStudyData(argv[1]);

  const auto start = Clock::now();

  auto faults   = study_data.model_data.bus_fault;
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
  std::cout << "\n\nComplete in " << dur << "\n";

  TestStatus status;
  for (std::size_t i = 0; i < stat_vec.size(); ++i)
  {
    status *= stat_vec[i];
    if (!stat_vec[i])
    {
      std::cout << "Study failed for fault: "
                << faults[i].disambiguation_string << '\n';
    }
  }

  return status.get();
}
