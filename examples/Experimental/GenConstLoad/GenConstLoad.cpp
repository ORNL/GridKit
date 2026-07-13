
#include <iomanip>
#include <iostream>

#include <IpIpoptApplication.hpp>
#include <IpSolveStatistics.hpp>

#include <GridKit/Model/PowerFlow/Bus/BusPQ.hpp>
#include <GridKit/Model/PowerFlow/Generator4Governor/Generator4Governor.hpp>
#include <GridKit/Model/PowerFlow/Load/Load.hpp>
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

  // Create a bus
  BaseBus<double, size_t>* bus = new BusPQ<double, size_t>(1.0, 0.0);

  // Attach a generator to the bus and signal ports
  ModelEvaluatorImpl<double, size_t>* gen = new Generator4Governor<double, size_t>(bus, 0.8, 0.3);

  // Attach load to the bus
  ModelEvaluatorImpl<double, size_t>* load = new Load<double, size_t>(bus, 0.8, 0.3);

  // Create system model
  SystemModel<double, size_t>* model = new SystemModel<double, size_t>();
  model->addBus(bus);
  model->addComponent(gen);
  model->addComponent(load);

  // Create numerical integrator and configure it for the generator model
  Ida<double, size_t>* idas = new Ida<double, size_t>(model);

  model->allocate();

  double t_init  = 0.0;
  double t_final = 15.0;

  // setup simulation
  idas->setMaxSteps(1000);
  idas->setBackwardMaxSteps(1000);
  idas->configureSimulation();
  idas->configureAdjoint();
  idas->getDefaultInitialCondition();
  idas->initializeSimulation(t_init);
  idas->configureQuadrature();
  idas->initializeQuadrature();

  idas->runSimulationQuadrature(0.1, 2);
  idas->saveInitialCondition();

  // Set integration time for dynamic constrained optimization
  idas->setIntegrationTime(t_init, t_final, 250);

  // Guess optimization parameter values
  double T2 = 0.15;
  double K  = 16.0;

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
  ipoptApp->Options()->SetIntegerValue("print_level", 5);

  // Create interface to Ipopt solver
  Ipopt::SmartPtr<Ipopt::TNLP> ipoptDynamicObjectiveInterface =
      new IpoptInterface::DynamicObjective<double, size_t>(idas);

  auto* param = model->param().getData();

  // Initialize problem
  param[0] = T2;
  param[1] = K;

  // Solve the problem
  status = ipoptApp->OptimizeTNLP(ipoptDynamicObjectiveInterface);

  if (status == Ipopt::Solve_Succeeded)
  {
    // Print result
    std::cout << "\nSucess: The problem solved in "
              << ipoptApp->Statistics()->IterationCount()
              << " iterations!\n";
    std::cout << "Optimal value: T2 = "
              << param[0]
              << ", K = "
              << param[1] << "\n";
    std::cout << "The final value of the objective function G(T2,K) = "
              << ipoptApp->Statistics()->FinalObjective() << "\n\n";
  }

  // Compare results of the two optimization methods
  int retval =
      isEqual(ipoptApp->Statistics()->FinalObjective(), 1239.0, 10 * tol) ? 0 : 1;

  if (retval != 0)
  {
    std::cout << "The two results differ beyond solver tolerance!\n";
  }

  delete idas;
  delete gen;
  delete load;
  delete bus;
  delete model;

  return 0;
}
