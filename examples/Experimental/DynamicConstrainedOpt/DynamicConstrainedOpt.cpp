#include <iomanip>
#include <iostream>
#include <memory>

#include <IpIpoptApplication.hpp>
#include <IpSolveStatistics.hpp>

#include <GridKit/Model/PowerFlow/Bus/BusSlack.hpp>
#include <GridKit/Model/PowerFlow/Generator2/Generator2.hpp>
#include <GridKit/Model/PowerFlow/SystemModel.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Solver/Optimization/DynamicConstraint.hpp>
#include <GridKit/Solver/Optimization/DynamicObjective.hpp>
#include <GridKit/Testing/Testing.hpp>

int main()
{
  using namespace GridKit;
  using namespace AnalysisManager::Sundials;
  using namespace AnalysisManager;
  using namespace GridKit::Testing;

  // Create an infinite bus
  BusSlack<double, size_t> bus(1.0, 0.0);

  // Attach a generator to that bus
  Generator2<double, size_t> gen(&bus);

  // Create a system model
  SystemModel<double, size_t> model;
  model.addBus(&bus);
  model.addComponent(&gen);

  // allocate model components
  model.allocate();

  // Create numerical integrator and configure it for the generator model
  Ida<double, size_t> idas(&model);

  double t_init  = 0.0;
  double t_final = 20.0;

  // setup simulation
  idas.configureSimulation();
  idas.setTolerance(1e-7, 1e-9);
  idas.configureAdjoint();
  idas.getDefaultInitialCondition();
  idas.initializeSimulation(t_init);
  idas.configureQuadrature();
  idas.initializeQuadrature();

  double t_fault = 0.02;
  double t_clear = 0.06;
  idas.runSimulation(t_fault);
  // create initial condition after a fault
  {
    gen.V() = 0.0;
    idas.runSimulation(t_clear, 2);
    gen.V()     = 1.0;
    gen.theta() = -0.01;
    idas.saveInitialCondition();
  }

  // Set integration time for dynamic constrained optimization
  idas.setIntegrationTime(t_init, t_final, 100);

  // Guess optimization parameter value
  double Pm = 0.7;

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

  // Set solver tolerance
  const double tol = 1e-4;

  // Configure Ipopt application
  ipoptApp->Options()->SetStringValue("hessian_approximation", "limited-memory");
  ipoptApp->Options()->SetNumericValue("tol", tol);
  ipoptApp->Options()->SetIntegerValue("print_level", 0);

  // Create dynamic objective interface to Ipopt solver
  Ipopt::SmartPtr<Ipopt::TNLP> ipoptDynamicObjectiveInterface =
      new IpoptInterface::DynamicObjective<double, size_t>(&idas);

  // Initialize problem
  model.param()[0] = Pm;

  // Solve the problem
  status = ipoptApp->OptimizeTNLP(ipoptDynamicObjectiveInterface);
  std::cout << "\n\nProblem formulated as dynamic objective optimization ...\n";

  if (status == Ipopt::Solve_Succeeded)
  {
    // Print result
    std::cout << "\nSucess:\n The problem solved in "
              << ipoptApp->Statistics()->IterationCount() << " iterations!\n"
              << " Optimal value of Pm = " << model.param()[0] << "\n"
              << " The final value of the objective function G(Pm) = "
              << ipoptApp->Statistics()->FinalObjective() << "\n\n";
  }

  // Store dynamic objective optimization results
  double* results = new double[model.sizeParams()];
  for (unsigned i = 0; i < model.sizeParams(); ++i)
  {
    results[i] = model.param()[i];
  }

  // Create dynamic constraint interface to Ipopt solver
  Ipopt::SmartPtr<Ipopt::TNLP> ipoptDynamicConstraintInterface =
      new IpoptInterface::DynamicConstraint<double, size_t>(&idas);

  // Initialize problem
  model.param()[0] = Pm;

  // Solve the problem
  status = ipoptApp->OptimizeTNLP(ipoptDynamicConstraintInterface);
  std::cout << "\n\nProblem formulated as dynamic constraint optimization ...\n";

  if (status == Ipopt::Solve_Succeeded)
  {
    // Print result
    std::cout << "\nSucess:\n The problem solved in "
              << ipoptApp->Statistics()->IterationCount() << " iterations!\n"
              << " Optimal value of Pm = " << model.param()[0] << "\n"
              << " The final value of the objective function G(Pm) = "
              << ipoptApp->Statistics()->FinalObjective() << "\n\n";
  }

  // Compare results of the two optimization methods
  int retval = 0;
  for (unsigned i = 0; i < model.sizeParams(); ++i)
  {
    if (!isEqual(results[i], model.param()[i], 10 * tol))
      --retval;
  }

  if (retval < 0)
  {
    std::cout << "The two results differ beyond solver tolerance!\n";
  }

  delete[] results;
  return retval;
}
