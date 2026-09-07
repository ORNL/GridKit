#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

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
  const double tol             = 1e-6;
  const double integration_tol = 1e-9;
  const double optimizer_tol   = 1e-8;
  const double objective_scale = 1e4;

  std::istringstream input_data(lookup_table);
  GridKit::setLookupTable(gen.getLookupTable(), input_data, t_init, t_final);

  std::cout << "Performing parameter estimation with respect to data\nfrom "
            << "t_init = " << t_init << " to t_final = " << t_final << "\n";

  // Resolve the squared-error objective and its gradient near the minimum.
  idas.setTolerance(integration_tol, integration_tol);
  idas.setBackwardTolerance(integration_tol, integration_tol);
  idas.setQuadratureTolerance(1e-8, 1e-13);
  idas.setBackwardQuadratureTolerance(1e-8, 1e-12);
  idas.setBackwardMaxSteps(10000);
  model.initialize();
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
    idas.runSimulation(t_clear, (t_clear - t_init) / 20.0);
    gen.V() = 1.0;
    idas.saveInitialCondition();
  }

  // Set monitoring interval for dynamic constrained optimization
  double dt_monitor = (t_final - t_init) / 100.0;

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

  // Scale the objective and specify unscaled stationarity and feasibility limits.
  ipoptApp->Options()->SetStringValue("hessian_approximation", "limited-memory");
  ipoptApp->Options()->SetNumericValue("tol", objective_scale * optimizer_tol);
  ipoptApp->Options()->SetStringValue("nlp_scaling_method", "none");
  ipoptApp->Options()->SetNumericValue("obj_scaling_factor", objective_scale);
  ipoptApp->Options()->SetNumericValue("dual_inf_tol", optimizer_tol);
  ipoptApp->Options()->SetNumericValue("constr_viol_tol", 1e-12);
  ipoptApp->Options()->SetIntegerValue("print_level", 0);

  // Create dynamic objective interface to Ipopt solver
  Ipopt::SmartPtr<Ipopt::TNLP> ipoptDynamicObjectiveInterface =
      new IpoptInterface::DynamicObjective<double, size_t>(&idas, t_init, t_final, dt_monitor);

  // Solve the problem
  status = ipoptApp->OptimizeTNLP(ipoptDynamicObjectiveInterface);
  if (status != Ipopt::Solve_Succeeded)
  {
    std::cerr << "Optimization failed with status " << status << "\n";
    return static_cast<int>(status);
  }
  std::cout << "\n\nProblem formulated as dynamic objective optimiztion ...\n";

  // Print result
  std::cout << "\nSucess:\n The problem solved in "
            << ipoptApp->Statistics()->IterationCount() << " iterations!\n"
            << " Optimal value of H = " << param[0] << "\n"
            << " The final value of the objective function G(H) = "
            << ipoptApp->Statistics()->FinalObjective() << "\n\n";

  // Store dynamic objective optimization results
  std::vector<double> results(model.sizeParams());
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
  if (status != Ipopt::Solve_Succeeded)
  {
    std::cerr << "Optimization failed with status " << status << "\n";
    return static_cast<int>(status);
  }
  std::cout << "\n\nProblem formulated as dynamic constraint optimiztion ...\n";

  // Print result
  std::cout << "\nSucess:\n The problem solved in "
            << ipoptApp->Statistics()->IterationCount() << " iterations!\n"
            << " Optimal value of H = " << param[0] << "\n"
            << " The final value of the objective function G(H) = "
            << ipoptApp->Statistics()->FinalObjective() << "\n\n";

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

  return retval;
}
