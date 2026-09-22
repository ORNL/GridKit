#include <chrono>
#include <exception>
#include <filesystem>
#include <fstream>
#include <future>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/PhasorDynamics/BusFault/BusFault.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

#include "AnalysisUtilities.hpp"
#include "SolverTrace.hpp"

using Clock = std::chrono::high_resolution_clock;
using Dur   = std::chrono::duration<double>;

using Log = GridKit::Utilities::Logger;

using namespace AnalysisManager::Sundials;
using namespace GridKit::PhasorDynamics;
using namespace GridKit::Testing;

using scalar_type = double;
using real_type   = double;
using index_type  = size_t;

struct StudyResult
{
  TestStatus            status{true};
  IdaStats              stats;
  std::vector<IdaStats> segments;
  std::string           diagnostic;
};

StudyResult runStudy(StudyData study_data)
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

  SolverTrace trace(study_data.solver_trace_file);
  const auto  accepted_step_callback = trace.callback(ida);

  StudyResult result;
  const auto  record_stats = [&]()
  {
    if (!study_data.contingency_stats_file.empty())
    {
      const auto stats  = ida.getStats();
      result.stats     += stats;
      result.segments.push_back(stats);
    }
  };

  using EventType = SystemEvent::Type;

  // Initilize simultation for first run
  real_type dt_monitor = study_data.dt_monitor;
  real_type final_time = study_data.tmax;
  ida.initializeSimulation(0.0, false);
  trace.record("init", 0.0, ida);

  for (const auto& event : study_data.events)
  {
    // Run to event time
    ida.runSimulation(event.time, dt_monitor, {}, accepted_step_callback);
    trace.finish(event.time, ida);
    record_stats();

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
    trace.record("init", event.time, ida);
  }

  // Run to final time
  ida.runSimulation(final_time, dt_monitor, {}, accepted_step_callback);
  trace.finish(final_time, ida);
  record_stats();

  // Stop the variable monitor
  sys.stopMonitor();

  result.status = checkErrors(study_data, false);
  if (result.status)
  {
    trace.write();
  }
  return result;
}

StudyResult singleFaultStudy(std::size_t fault_id, StudyData study_data)
{
  // Change id in schedule to current fault id
  for (auto& event : study_data.events)
  {
    event.element_id = fault_id;
  }

  // Make distinct output files
  for (auto& sink : study_data.model_data.monitor_sink)
  {
    auto path      = std::filesystem::path(sink.file_name);
    auto ext       = path.extension().string();
    auto name      = path.stem().string();
    sink.file_name = name + "_" + std::to_string(fault_id) + ext;
  }

  if (!study_data.solver_trace_file.empty() && study_data.model_data.bus_fault.size() > 1)
  {
    const auto path              = study_data.solver_trace_file;
    study_data.solver_trace_file = path.parent_path() / (path.stem().string() + "_" + std::to_string(fault_id) + path.extension().string());
  }

  try
  {
    return runStudy(study_data);
  }
  catch (const std::exception& error)
  {
    Log::warning() << "exception caught at fault id: " << fault_id << ": " << error.what() << std::endl;
    return {{false}, {}, {}, error.what()};
  }
  catch (...)
  {
    Log::warning() << "exception caught at fault id: " << fault_id << std::endl;
    return {{false}, {}, {}, "Unknown exception"};
  }
}

void runStudySerial(const StudyData& study_data, std::vector<StudyResult>& stat_vec)
{
  for (std::size_t i = 0; i < study_data.model_data.bus_fault.size(); ++i)
  {
    auto stat   = singleFaultStudy(i, study_data);
    stat_vec[i] = stat;
  }
}

#ifdef GRIDKIT_ENABLE_THREADS
void runStudyAsync(const StudyData& study_data, std::vector<StudyResult>& stat_vec)
{
  auto n_faults = study_data.model_data.bus_fault.size();

  std::vector<std::future<StudyResult>> futures;
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
void runStudyOpenMP(const StudyData& study_data, std::vector<StudyResult>& stat_vec)
{
  auto n_faults = study_data.model_data.bus_fault.size();
#pragma omp parallel for
  for (std::size_t i = 0; i < n_faults; ++i)
  {
    auto stat   = singleFaultStudy(i, study_data);
    stat_vec[i] = stat;
  }
}
#endif

json statisticsJson(const IdaStats& stats)
{
  return {
      {"accepted_steps", stats.num_steps_},
      {"jacobian_evals", stats.num_jacobian_evals_},
      {"residual_evals", stats.num_residual_evals_},
      {"linear_solver_setups", stats.num_linear_decompositions_},
      {"error_test_failures", stats.num_error_test_fails_},
      {"nonlinear_iterations", stats.num_nonlinear_iters_},
      {"nonlinear_convergence_failures", stats.num_nonlinear_convergence_fails_}};
}

void writeStatistics(const StudyData& study, const std::vector<StudyResult>& results)
{
  json records = json::array();
  for (std::size_t i = 0; i < results.size(); ++i)
  {
    const auto& result = results[i];
    const auto& fault  = study.model_data.bus_fault[i];
    json        record = {
        {"fault_id", i},
        {"bus", fault.buses.at(BusFaultBuses::bus)},
        {"status", result.status ? "ok" : "failed"},
        {"diagnostic", result.diagnostic},
        {"stats", nullptr},
        {"segments", json::array()}};
    if (result.status)
    {
      record["stats"] = statisticsJson(result.stats);
      double start    = 0.0;
      for (std::size_t segment = 0; segment < result.segments.size(); ++segment)
      {
        const double end   = segment < study.events.size() ? study.events[segment].time : study.tmax;
        auto         item  = statisticsJson(result.segments[segment]);
        item["start_time"] = start;
        item["end_time"]   = end;
        record["segments"].push_back(item);
        start = end;
      }
    }
    records.push_back(record);
  }

  const json document = {
      {"schema_version", 1},
      {"statistics_scope", "IDA segment totals, including event consistent-initialization work"},
      {"initial_consistent_ic", false},
      {"records", records}};
  std::ofstream output;
  output.exceptions(std::ios::failbit | std::ios::badbit);
  output.open(study.contingency_stats_file);
  output << document.dump(2) << '\n';
  output.close();
}

int runApplication(int argc, const char* argv[])
{
  // Study file
  checkCommandLine(argc, "ContingencyAnalysis");
  auto study_data = parseStudyData(argv[1]);

  GridKit::Math::MU<real_type> = study_data.mu;

  const auto start = Clock::now();

  auto faults   = study_data.model_data.bus_fault;
  auto stat_vec = std::vector<StudyResult>(faults.size());

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
    status *= stat_vec[i].status;
    if (!stat_vec[i].status)
    {
      std::cout << "Study failed for fault: "
                << faults[i].disambiguation_string << '\n';
    }
  }

  if (!study_data.contingency_stats_file.empty())
  {
    writeStatistics(study_data, stat_vec);
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
