#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

#include <IpIpoptApplication.hpp>
#include <IpSolveStatistics.hpp>

#include <GridKit/Model/PowerFlow/Bus/BusSlack.hpp>
#include <GridKit/Model/PowerFlow/Generator4Param/Generator4Param.hpp>
#include <GridKit/Model/PowerFlow/SystemModel.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Solver/Optimization/DynamicConstraint.hpp>
#include <GridKit/Solver/Optimization/DynamicObjective.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/FileIO.hpp>

#include "lookup_table.hpp"

int main()
{
  using namespace GridKit;
  using namespace AnalysisManager::Sundials;
  using namespace AnalysisManager;
  using namespace GridKit::Testing;

  // Create an infinite bus
  BusSlack<double, size_t> bus(1.0, 0.0);

  // Attach a generator to that bus
  Generator4Param<double, size_t> gen(&bus);

  // Create a system model
  SystemModel<double, size_t> model;
  model.addBus(&bus);
  model.addComponent(&gen);

  // allocate model components
  model.allocate();

  // Create numerical integrator and configure it for the generator model
  Ida<double, size_t> idas(&model);

  double t_init  = -1.0;
  double t_final = -1.0;

  // Set solver tolerance
  const double tol = 1e-6;

  std::istringstream input_data(lookup_table);
  GridKit::setLookupTable(gen.getLookupTable(), input_data, t_init, t_final);

  std::cout << "Performing parameter estimation with respect to data\nfrom "
            << "t_init = " << t_init << " to t_final = " << t_final << "\n";

  // setup simulation
  idas.setTolerance(0.1 * tol);
  idas.setBackwardTolerance(0.1 * tol);
  idas.setQuadratureTolerance(0.1 * tol);
  idas.setBackwardQuadratureTolerance(0.1 * tol);
  idas.configureSimulation();
  idas.configureAdjoint();
  idas.getDefaultInitialCondition();
  idas.initializeSimulation(t_init);
  idas.configureQuadrature();
  idas.initializeQuadrature();

  double t_fault = 0.1;
  double t_clear = 0.1;
  idas.runSimulation(t_fault);
  idas.saveInitialCondition();
  // create initial condition after a fault
  {
    idas.getSavedInitialCondition();
    gen.V() = 0.0;
    idas.initializeSimulation(t_init);
    idas.runSimulation(t_clear, 20);
    gen.V() = 1.0;
    idas.saveInitialCondition();
  }

  // Set monitoring interval for dynamic constrained optimization
  double dt_monitor = (t_final - t_init) / 100.0;

  auto* param = model.param().getData();

  auto* param = model.param().getData();

  // Guess value of inertia coefficient
  param[0] = 3.0;
  model.param().setDataUpdated();

  // Create an instance of the IpoptApplication
  Ipopt::SmartPtr<Ipopt::IpoptApplication> ipoptApp = IpoptApplicationFactory();

  // Initialize the IpoptApplication and process the options
  Ipopt::ApplicationReturnStatus status;
  status = ipoptApp->Initialize();
  if (status != Ipopt::Solve_Succeeded)
  {
    std::cout << "\n\n*** Initialization failed! ***\n\n";
    return (int) status;
  }

  // Configure Ipopt application
  ipoptApp->Options()->SetStringValue("hessian_approximation", "limited-memory");
  ipoptApp->Options()->SetNumericValue("tol", tol);
  ipoptApp->Options()->SetIntegerValue("print_level", 0);

  // Create dynamic objective interface to Ipopt solver
  Ipopt::SmartPtr<Ipopt::TNLP> ipoptDynamicObjectiveInterface =
      new IpoptInterface::DynamicObjective<double, size_t>(&idas, t_init, t_final, dt_monitor);

  // Solve the problem
  status = ipoptApp->OptimizeTNLP(ipoptDynamicObjectiveInterface);
  std::cout << "\n\nProblem formulated as dynamic objective optimiztion ...\n";

  if (status == Ipopt::Solve_Succeeded)
  {
    // Print result
    std::cout << "\nSucess:\n The problem solved in "
              << ipoptApp->Statistics()->IterationCount() << " iterations!\n"
              << " Optimal value of H = " << param[0] << "\n"
              << " The final value of the objective function G(H) = "
              << ipoptApp->Statistics()->FinalObjective() << "\n\n";
  }

  // Store dynamic objective optimization results
  double* results = new double[model.sizeParams()];
  for (unsigned i = 0; i < model.sizeParams(); ++i)
  {
    results[i] = param[i];
  }

  // Guess value of inertia coefficient
  param[0] = 3.0;

  // Create dynamic constraint interface to Ipopt solver
  Ipopt::SmartPtr<Ipopt::TNLP> ipoptDynamicConstraintInterface =
      new IpoptInterface::DynamicConstraint<double, size_t>(&idas, t_init, t_final, dt_monitor);

  // Solve the problem
  status = ipoptApp->OptimizeTNLP(ipoptDynamicConstraintInterface);
  std::cout << "\n\nProblem formulated as dynamic constraint optimiztion ...\n";

  if (status == Ipopt::Solve_Succeeded)
  {
    // Print result
    std::cout << "\nSucess:\n The problem solved in "
              << ipoptApp->Statistics()->IterationCount() << " iterations!\n"
              << " Optimal value of H = " << param[0] << "\n"
              << " The final value of the objective function G(H) = "
              << ipoptApp->Statistics()->FinalObjective() << "\n\n";
  }

  // Compare results of the two optimization methods
  int retval = 0;
  for (unsigned i = 0; i < model.sizeParams(); ++i)
  {
    if (!isEqual(results[i], param[i], 100 * tol))
      --retval;
  }

  if (retval < 0)
  {
    std::cout << "The two results differ beyond solver tolerance!\n";
  }

  delete[] results;
  return retval;
}
