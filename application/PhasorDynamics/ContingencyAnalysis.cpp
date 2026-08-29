#include <chrono>
#include <filesystem>
#include <future>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <GridKit/Model/PhasorDynamics/BusFault/BusFault.hpp>
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
  ida.setOptions(study_data.ida);
  ida.setConsistentICType(study_data.consistent_ic_type);
  ida.configureSimulation();

  using EventType = SystemEvent::Type;

  // Initilize simultation for first run
  real_type dt_monitor = study_data.dt_monitor;
  real_type final_time = study_data.tmax;
  ida.initializeSimulation(0.0, false);

  for (const auto& event : study_data.events)
  {
    // Run to event time
    ida.runSimulation(event.time, dt_monitor);

    // Set up run for event (to start at event time)
    switch (event.type)
    {
    case EventType::FAULT_ON:
      sys.getBusFault(event.element_id)->setStatus(true);
      break;
    case EventType::FAULT_OFF:
      sys.getBusFault(event.element_id)->setStatus(false);
      break;
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

TestStatus singleFaultStudy(std::size_t bus_id, StudyData study_data)
{
  // Relocate the single fault template and direct every fault event to it.
  study_data.model_data.bus_fault.front().buses[BusFaultBuses::bus] = bus_id;
  for (auto& event : study_data.events)
  {
    event.element_id = 0;
  }

  // Make distinct output files
  for (auto& sink : study_data.model_data.monitor_sink)
  {
    auto path      = std::filesystem::path(sink.file_name);
    auto ext       = path.extension().string();
    auto name      = path.stem().string();
    sink.file_name = name + "_bus" + std::to_string(bus_id) + ext;
  }

  try
  {
    return runStudy(study_data);
  }
  catch (...)
  {
    Log::warning() << "exception caught at bus: " << bus_id << std::endl;
    return {false};
  }
}

void runStudySerial(const StudyData&               study_data,
                    const std::vector<index_type>& bus_ids,
                    std::vector<TestStatus>&       stat_vec)
{
  for (std::size_t i = 0; i < bus_ids.size(); ++i)
  {
    auto stat   = singleFaultStudy(bus_ids[i], study_data);
    stat_vec[i] = stat;
  }
}

#ifdef GRIDKIT_ENABLE_THREADS
void runStudyAsync(const StudyData&               study_data,
                   const std::vector<index_type>& bus_ids,
                   std::vector<TestStatus>&       stat_vec)
{
  std::vector<std::future<TestStatus>> futures;
  futures.reserve(bus_ids.size());
  for (const auto bus_id : bus_ids)
  {
    futures.emplace_back(
        std::async(std::launch::async, singleFaultStudy, bus_id, study_data));
  }

  for (std::size_t i = 0; i < bus_ids.size(); ++i)
  {
    auto stat   = futures[i].get();
    stat_vec[i] = stat;
  }
}
#endif

#ifdef _OPENMP
void runStudyOpenMP(const StudyData&               study_data,
                    const std::vector<index_type>& bus_ids,
                    std::vector<TestStatus>&       stat_vec)
{
#pragma omp parallel for
  for (std::size_t i = 0; i < bus_ids.size(); ++i)
  {
    auto stat   = singleFaultStudy(bus_ids[i], study_data);
    stat_vec[i] = stat;
  }
}
#endif

int main(int argc, const char* argv[])
{
  // Study file
  checkCommandLine(argc, "ContingencyAnalysis");
  auto study_data = parseStudyData(argv[1]);

  if (study_data.model_data.bus_fault.empty())
  {
    Log::error() << "ContingencyAnalysis requires a BusFault in the system model\n";
    return 1;
  }

  std::vector<index_type> bus_ids;
  bus_ids.reserve(study_data.model_data.bus.size());
  for (const auto& bus : study_data.model_data.bus)
  {
    bus_ids.push_back(bus.bus_id);
  }

  const auto start = Clock::now();

  auto stat_vec = std::vector<TestStatus>(bus_ids.size(), true);

  // Use std::async if threads are available (so far, std::async has out-performed OpenMP)
  // Otherwise, use OpenMP if available
  // Fall back to serial if neither threads or OpenMP are available
#if defined(GRIDKIT_ENABLE_THREADS)
  runStudyAsync(study_data, bus_ids, stat_vec);
#elif defined(_OPENMP)
  runStudyOpenMP(study_data, bus_ids, stat_vec);
#else
  runStudySerial(study_data, bus_ids, stat_vec);
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
      std::cout << "Study failed for bus: " << bus_ids[i] << '\n';
    }
  }

  return status.get();
}
