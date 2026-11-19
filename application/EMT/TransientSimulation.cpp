/**
 * @file TransientSimulation.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief EMT transient simulation of a study file.
 */

#include <ctime>
#include <exception>

#include <GridKit/Model/EMT/StateJSONParser.hpp>
#include <GridKit/Model/EMT/SystemModel.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>

#include "AnalysisUtilities.hpp"

using namespace AnalysisManager::Sundials;
using namespace GridKit::EMT;

using scalar_type = double;
using real_type   = double;
using index_type  = size_t;

int runApplication(int argc, const char* argv[])
{
  Log::raiseVerbosity(Log::SUMMARY);

  checkCommandLine(argc, "TransientSimulation");
  auto study = parseStudyData(argv[1]);

  SystemModel<scalar_type, index_type> sys(study.model_data);
  if (sys.allocate() != 0)
  {
    return 1;
  }
  applyState(study.state, sys);

  Ida<scalar_type, index_type> ida(&sys);
  ida.setTolerance(study.rel_tol, study.abs_tol);
  ida.setFixedStep(study.dt_fixed);
  ida.setMaxSteps(study.max_steps);
  ida.setMaxOrder(study.max_order);
  ida.setConsistentICType(study.consistent_ic_type);
  if (ida.configureSimulation() != 0 || ida.initializeSimulation(0.0) != 0)
  {
    return 1;
  }

  real_type start = static_cast<real_type>(clock());

  // Events at one time take effect together, with one restart
  auto event = study.events.cbegin();
  while (event != study.events.cend())
  {
    const auto time = event->at("time").get<real_type>();
    ida.runSimulation(time, study.dt_monitor);
    do
    {
      applyState(*event, sys);
      ++event;
    } while (event != study.events.cend() && event->at("time").get<real_type>() == time);

    // The solver restarts from the state the events have set
    if (ida.getDefaultInitialCondition() != 0 || ida.initializeSimulation(time) != 0)
    {
      return 1;
    }
  }
  ida.runSimulation(study.tmax, study.dt_monitor);

  real_type stop = static_cast<real_type>(clock());

  sys.stopMonitor();

  Log::summary() << "Complete in " << (stop - start) / CLOCKS_PER_SEC << " seconds\n";
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
    Log::error() << "TransientSimulation failed: " << error.what() << '\n';
  }

  return 1;
}
