#include <chrono>
#include <exception>
#include <filesystem>
#include <future>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <GridKit/Model/PhasorDynamics/Branch/Branch.hpp>
#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

#include "CoSimulation.hpp"

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
  return runStandaloneStudy(study_data, false, false);
}

TestStatus singleFaultStudy(std::size_t fault_id, StudyData study_data)
{
  // Retarget fault events while preserving branch-switching events
  for (auto& event : study_data.events)
  {
    if (event.type == SystemEvent::Type::FAULT_ON
        || event.type == SystemEvent::Type::FAULT_OFF)
    {
      event.bus = study_data.model_data.bus[fault_id].bus_id;
    }
  }

  // Make distinct output files
  suffixOutputs(study_data, "_" + std::to_string(fault_id));

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

void runStudySerial(
  const StudyData& study_data, std::vector<TestStatus>& stat_vec
)
{
  for (std::size_t i = 0; i < study_data.model_data.bus.size(); ++i)
  {
    auto stat   = singleFaultStudy(i, study_data);
    stat_vec[i] = stat;
  }
}

#ifdef GRIDKIT_ENABLE_THREADS
void runStudyAsync(
  const StudyData& study_data, std::vector<TestStatus>& stat_vec
)
{
  auto n_faults = study_data.model_data.bus.size();

  std::vector<std::future<TestStatus>> futures;
  futures.reserve(n_faults);
  for (std::size_t i = 0; i < n_faults; ++i)
  {
    futures.emplace_back(
      std::async(std::launch::async, singleFaultStudy, i, study_data)
    );
  }

  for (std::size_t i = 0; i < n_faults; ++i)
  {
    auto stat   = futures[i].get();
    stat_vec[i] = stat;
  }
}
#endif

#ifdef _OPENMP
void runStudyOpenMP(
  const StudyData& study_data, std::vector<TestStatus>& stat_vec
)
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
  const auto options    = parseSimulationOptions(argc, argv, true);
  auto       study_data = parseStudyData(options.study_file);

  const auto start = Clock::now();

  auto faults   = study_data.model_data.bus;
  auto stat_vec = std::vector<TestStatus>(faults.size(), true);

#ifdef GRIDKIT_ENABLE_ZMQ
  if (options.coSimulated())
  {
    std::size_t next_fault = 0;
    std::size_t completed  = 0;
    const auto  next       = [&]() -> std::optional<StudyData>
    {
      if (next_fault == faults.size())
      {
        return std::nullopt;
      }

      const auto i    = next_fault++;
      auto       data = study_data;
      for (auto& event : data.events)
      {
        if (event.type == SystemEvent::Type::FAULT_ON
            || event.type == SystemEvent::Type::FAULT_OFF)
        {
          event.bus = faults[i].bus_id;
        }
      }

      suffixOutputs(data, "_" + std::to_string(i));
      return data;
    };
    if (runCoupledClientStudies(
          options,
          next,
          [&](TestStatus result)
          { stat_vec[completed++] = result; },
          true
        )
        == GridKit::CoSim::RunResult::STOPPED)
    {
      return 1;
    }
  }
  else
#endif
  {
    // Use std::async if threads are available
    // Otherwise, use OpenMP if available
    // Fall back to serial if neither threads or OpenMP are available
#if defined(GRIDKIT_ENABLE_THREADS)
    runStudyAsync(study_data, stat_vec);
#elif defined(_OPENMP)
    runStudyOpenMP(study_data, stat_vec);
#else
    runStudySerial(study_data, stat_vec);
#endif
  }

  const auto stop = Clock::now();
  const auto dur  = std::chrono::duration<double>(stop - start);
  Log::summary() << "Complete in " << dur << "\n";

  TestStatus status;
  for (std::size_t i = 0; i < stat_vec.size(); ++i)
  {
    status *= stat_vec[i];
    if (!stat_vec[i])
    {
      Log::error() << "Study failed for fault: " << faults[i].name << '\n';
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
