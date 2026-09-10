#include <chrono>
#include <exception>

#include "CoSimulation.hpp"

using namespace GridKit::PhasorDynamics;

int runApplication(int argc, const char* argv[])
{
  Log::raiseVerbosity(Log::SUMMARY);
  const auto                   options = parseSimulationOptions(argc, argv);
  const auto                   study   = parseStudyData(options.study_file);
  const auto                   start   = std::chrono::steady_clock::now();
  GridKit::Testing::TestStatus status;
#ifdef GRIDKIT_ENABLE_ZMQ
  if (options.coSimulated())
  {
    status = runCoupledStudy(study, options);
  }
  else
#endif
  {
    status = runStandaloneStudy(study);
  }

  const auto elapsed =
    std::chrono::duration<double>(std::chrono::steady_clock::now() - start);
  Log::summary() << "Complete in " << elapsed.count() << " seconds\n";
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
