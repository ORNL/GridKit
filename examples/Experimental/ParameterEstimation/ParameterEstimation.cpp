#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

#include "lookup_table.hpp"
#include <IpIpoptApplication.hpp>
#include <IpSolveStatistics.hpp>
#include <Model/PowerFlow/Bus/BusSlack.hpp>
#include <Model/PowerFlow/Generator4Param/Generator4Param.hpp>
#include <Model/PowerFlow/SystemModel.hpp>
#include <Solver/Dynamic/Ida.hpp>
#include <Solver/Optimization/DynamicConstraint.hpp>
#include <Solver/Optimization/DynamicObjective.hpp>
#include <Utilities/FileIO.hpp>
#include <Utilities/Testing.hpp>

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

  std::istringstream input_data(lookup_table);
  GridKit::setLookupTable(gen.getLookupTable(), input_data, t_init, t_final);

  std::cout << "Performing parameter estimation with respect to data\nfrom "
            << "t_init = " << t_init << " to t_final = " << t_final << "\n";

  // setup simulation
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
    idas.initializeSimulation(t_init);
    gen.V() = 0.0;
    idas.runSimulation(t_clear, 20);
    gen.V() = 1.0;
    idas.saveInitialCondition();
  }

  // Set integration time for dynamic constrained optimization
  idas.setIntegrationTime(t_init, t_final, 100);

  // Guess value of inertia coefficient
  model.param()[0] = 3.0;

  // Create an instance of the IpoptApplication
  Ipopt::SmartPtr<Ipopt::IpoptApplication> ipoptApp = IpoptApplicationFactory();

  // Set solver tolerance
  const double tol = 1e-5;

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
      new IpoptInterface::DynamicObjective<double, size_t>(&idas);

  // Solve the problem
  status = ipoptApp->OptimizeTNLP(ipoptDynamicObjectiveInterface);
  std::cout << "\n\nProblem formulated as dynamic objective optimiztion ...\n";

  if (status == Ipopt::Solve_Succeeded)
  {
    // Print result
    std::cout << "\nSucess:\n The problem solved in "
              << ipoptApp->Statistics()->IterationCount() << " iterations!\n"
              << " Optimal value of H = " << model.param()[0] << "\n"
              << " The final value of the objective function G(H) = "
              << ipoptApp->Statistics()->FinalObjective() << "\n\n";
  }

  // Store dynamic objective optimization results
  double* results = new double[model.sizeParams()];
  for (unsigned i = 0; i < model.sizeParams(); ++i)
  {
    results[i] = model.param()[i];
  }

  // Guess value of inertia coefficient
  model.param()[0] = 3.0;

  // Create dynamic constraint interface to Ipopt solver
  Ipopt::SmartPtr<Ipopt::TNLP> ipoptDynamicConstraintInterface =
      new IpoptInterface::DynamicConstraint<double, size_t>(&idas);

  // Solve the problem
  status = ipoptApp->OptimizeTNLP(ipoptDynamicConstraintInterface);
  std::cout << "\n\nProblem formulated as dynamic constraint optimiztion ...\n";

  if (status == Ipopt::Solve_Succeeded)
  {
    // Print result
    std::cout << "\nSucess:\n The problem solved in "
              << ipoptApp->Statistics()->IterationCount() << " iterations!\n"
              << " Optimal value of H = " << model.param()[0] << "\n"
              << " The final value of the objective function G(H) = "
              << ipoptApp->Statistics()->FinalObjective() << "\n\n";
  }

  // Compare results of the two optimization methods
  int retval = 0;
  for (unsigned i = 0; i < model.sizeParams(); ++i)
  {
    if (!isEqual(results[i], model.param()[i], 100 * tol))
      --retval;
  }

  if (retval < 0)
  {
    std::cout << "The two results differ beyond solver tolerance!\n";
  }

  delete[] results;
  return retval;
}
