#include <filesystem>
#include <fstream>

#include <GridKit/Model/EMT/SystemModel.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

#include "AnalysisUtilities.hpp"
#include "EventSchedule.hpp"
#include "StateMonitor.hpp"

using namespace AnalysisManager::Sundials;
using namespace GridKit::EMT;
using namespace GridKit::Testing;

using scalar_type = double;
using real_type   = double;
using index_type  = size_t;

int main(int argc, const char* argv[])
{
  // Study file
  checkCommandLine(argc, "EMTDynamicSimulation");
  auto study = parseStudyData(argv[1]);
  configureCommonMath<real_type>(study);

#ifndef GRIDKIT_ENABLE_SUNDIALS_SPARSE
  throw std::runtime_error("EMTDynamicSimulation requires SUNDIALS with sparse KLU support");
#endif

  // Instantiate system
  SystemModel<scalar_type, index_type>   sys(study.model_data);
  EventSchedule<scalar_type, index_type> events(sys, study);
  if (!sys.hasJacobian())
  {
    throw std::runtime_error("EMTDynamicSimulation requires a sparse model Jacobian; enable Enzyme");
  }
  sys.allocate();
  if (sys.initialize(study.state) != 0)
    throw std::runtime_error("EMT model initialization failed");

  // Set up simulation
  Ida<scalar_type, index_type> ida(&sys);
  ida.setTolerance(study.rel_tol, study.abs_tol);
  ida.setFixedStep(study.dt_fixed);
  ida.setMaxSteps(study.max_steps);
  ida.setConsistentICType(study.consistent_ic_type);
  events.configure(ida);
  std::cout << "Linear solver: SUNDIALS KLU (sparse)\n"
            << "DAE variables: " << sys.size()
            << ", Jacobian nonzeros: " << sys.getCsrJacobian()->getNnz()
            << ", mu: " << study.mu << std::endl;
  StateMonitor<scalar_type, index_type> state_monitor(sys, study);
  auto                                  record_state = [&](real_type time)
  { state_monitor.write(time); };

  // Start timer
  real_type start = static_cast<real_type>(clock());

  const auto total_stats = events.run(ida, record_state);

  real_type stop = static_cast<real_type>(clock());

  // Stop the variable monitor
  sys.stopMonitor();

  // Generate aggregate errors comparing variable output to reference solution
  TestStatus status = checkErrors(study);

  // Report run time
  std::cout << "\n\nComplete in " << (stop - start) / CLOCKS_PER_SEC << " seconds\n";
  std::cout << "IDA statistics: steps=" << total_stats.num_steps_
            << ", residual_evals=" << total_stats.num_residual_evals_
            << ", linear_setups=" << total_stats.num_linear_decompositions_
            << ", error_test_fails=" << total_stats.num_error_test_fails_
            << ", nonlinear_iters=" << total_stats.num_nonlinear_iters_
            << ", nonlinear_convergence_fails=" << total_stats.num_nonlinear_convergence_fails_
            << '\n';

  return status.get();
}
