#include <ctime>
#include <exception>
#include <iomanip>
#include <vector>

#include <GridKit/Model/PhasorDynamics/Branch/Branch.hpp>
#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/Partition/PartitionData.hpp>
#include <GridKit/Model/PhasorDynamics/Partition/PartitionPlacement.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Solver/Dynamic/Partitioned.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

#include "AnalysisUtilities.hpp"

using namespace AnalysisManager::Sundials;
using namespace GridKit::PhasorDynamics;
using namespace GridKit::Testing;

using scalar_type = double;
using real_type   = double;
using index_type  = size_t;

/**
 * @brief Runs a study with the partitioned DAE integrator.
 *
 * @todo Share the study loop with DynamicSimulation once that loop is generic
 * over the solver. Merge the two executables once Ida (IDAS) and Partitioned
 * (IDA) can be linked together.
 */
int runApplication(int argc, const char* argv[])
{
  // Print summaries, such as the run time, without lowering a higher verbosity
  Log::raiseVerbosity(Log::SUMMARY);

  // Study file
  checkCommandLine(argc, "PartitionedDynamicSimulation");
  auto study = parseStudyData(argv[1]);

  using EventType = SystemEvent::Type;

  int ret = 0;
  if (study.dt_fixed <= 0.0)
  {
    Log::error() << "The partitioned integrator requires a positive dt_fixed\n";
    ret = 1;
  }

  if (study.partition_file.empty())
  {
    Log::error() << "The partitioned integrator requires a partition_file\n";
    ret = 1;
  }

  for (size_t event_index = 0; event_index < study.events.size() && ret == 0; event_index++)
  {
    const auto& event = study.events[event_index];
    if ((event.type == EventType::FAULT_ON || event.type == EventType::FAULT_OFF)
        && event.bus == GridKit::INVALID_INDEX<index_type>)
    {
      Log::error() << "Fault event at time " << event.time
                   << " requires a bus number\n";
      ret = 1;
    }
  }

  int exit_code = 1;
  if (ret == 0)
  {
    // Split the case into partitions by the partition file
    const PartitionData<index_type> partition = parsePartitionData(study.partition_file);
    study.model_data =
        PartitionPlacement<real_type, index_type>::apply(study.model_data, partition);

    // Instantiate system
    SystemModel<scalar_type, index_type> sys(study.model_data);
    sys.allocate();

    // Each top-level subsystem is a partition and the remaining variables couple them
    const index_type  system_size = sys.size();
    Mask              coupling_mask(system_size, true);
    std::vector<Mask> component_masks(sys.numSubsystems(), Mask(system_size, false));
    for (index_type subsystem_index = 0; subsystem_index < sys.numSubsystems(); subsystem_index++)
    {
      auto*            subsystem = sys.getSubsystem(subsystem_index);
      const index_type first     = subsystem->offset();
      const index_type last      = first + subsystem->size();
      for (index_type variable = first; variable < last; variable++)
      {
        component_masks[subsystem_index][variable] = true;
        coupling_mask[variable]                    = false;
      }
    }

    Log::warning() << "The partitioned integrator ignores max_order and consistent_ic_type\n";
    if (!study.events.empty())
    {
      Log::warning() << "The partitioned integrator restarts at events"
                     << " without recomputing consistent algebraic variables\n";
    }

    // Set up simulation
    Partitioned<scalar_type, index_type> solver(&sys, component_masks, coupling_mask);
    solver.setTolerance(study.rel_tol, study.abs_tol);
    solver.setFixedStep(study.dt_fixed);
    solver.setMaxSteps(study.max_steps);
    solver.configureSimulation();

    // Start timer
    real_type start = static_cast<real_type>(clock());

    // Initialize simulation for first run
    auto      dt_monitor = study.dt_monitor;
    real_type final_time = study.tmax;
    solver.initializeSimulation(0.0);
    for (size_t event_index = 0; event_index < study.events.size() && ret == 0; event_index++)
    {
      const auto& event = study.events[event_index];

      // Run to event time
      solver.runSimulation(event.time, dt_monitor);

      // Set up run for event (to start at event time)
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
      }

      if (ret != 0)
      {
        Log::error() << "Fault event failed for bus " << event.bus
                     << " at time " << event.time << '\n';
      }
      else
      {
        // Re-initialize simulation at event time
        solver.initializeSimulation(event.time);
      }
    }

    if (ret == 0)
    {
      // Run to final time
      solver.runSimulation(final_time, dt_monitor);

      real_type stop = static_cast<real_type>(clock());

      // Stop the variable monitor
      sys.stopMonitor();

      // Generate aggregate errors comparing variable output to reference solution
      TestStatus status = checkErrors(study);

      // Report run time
      Log::summary() << "Complete in " << (stop - start) / CLOCKS_PER_SEC << " seconds\n";

      // Report the outcome of the integration
      PartitionedStats stats = solver.getStats();
      Log::summary() << "Partitioned integrator reached t = " << final_time
                     << " with " << sys.numSubsystems() << " partitions\n";
      int label_width = 30;
      int stat_width  = 12;
      Log::summary() << std::setw(label_width) << "Outer steps" << " : "
                     << std::setw(stat_width) << stats.num_steps_ << '\n';
      Log::summary() << std::setw(label_width) << "Coupling nonlinear iterations" << " : "
                     << std::setw(stat_width) << stats.num_coupling_nonlinear_iters_ << '\n';
      Log::summary() << std::setw(label_width) << "Coupling nonlinear failures" << " : "
                     << std::setw(stat_width) << stats.num_coupling_nonlinear_failures_ << '\n';
      if (!study.events.empty())
      {
        Log::summary() << "Each event restarts the integrator,"
                       << " so these counts cover only the run after the last event\n";
      }
      if (study.output_file.empty() || study.reference_file.empty())
      {
        Log::summary() << "The solution was not compared to a reference\n";
      }

      exit_code = status.get();
    }
  }

  return exit_code;
}

int main(int argc, const char* argv[])
{
  int exit_code = 1;
  try
  {
    exit_code = runApplication(argc, argv);
  }
  catch (const std::exception& error)
  {
    Log::error() << "PartitionedDynamicSimulation failed: " << error.what() << '\n';
  }

  return exit_code;
}
