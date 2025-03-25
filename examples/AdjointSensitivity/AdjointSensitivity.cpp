
#include <iomanip>
#include <iostream>

#include <Model/PowerFlow/Bus/BusSlack.hpp>
#include <Model/PowerFlow/Generator4/Generator4.hpp>
#include <Model/PowerFlow/SystemModel.hpp>
#include <Solver/Dynamic/Ida.hpp>
#include <Utilities/Testing.hpp>

/*
 * Compute gradient of an objective function expressed as an integral over
 * system trajectory. The gradient is computed numerically and using
 * adjoint sensitivity analysis.
 *
 * The test case is a 4th order generator connected to an infinite bus.
 * The objective function is total frequency deviation computed over
 * system trajectory after generator short circuit fault.
 *
 */
int main()
{
  using namespace GridKit;
  using namespace AnalysisManager::Sundials;
  using namespace AnalysisManager;
  using namespace GridKit::Testing;

  // Create an infinite bus
  BaseBus<double, size_t>* bus = new BusSlack<double, size_t>(1.0, 0.0);

  // Attach a generator to that bus
  Generator4<double, size_t>* gen = new Generator4<double, size_t>(bus, 0.8, 0.3);

  // Create a system model
  SystemModel<double, size_t>* model = new SystemModel<double, size_t>();
  model->addBus(bus);
  model->addComponent(gen);

  // allocate model components
  model->allocate();

  // Create numerical integrator and configure it for the generator model
  Ida<double, size_t>* idas = new Ida<double, size_t>(model);

  double t_init  = 0.0;
  double t_final = 15.0;

  // setup simulation
  idas->configureSimulation();
  idas->configureAdjoint();
  idas->getDefaultInitialCondition();
  idas->initializeSimulation(t_init);
  idas->configureQuadrature();
  idas->initializeQuadrature();

  idas->runSimulation(0.1, 2);
  idas->saveInitialCondition();

  // create initial condition after a fault
  {
    idas->getSavedInitialCondition();
    idas->initializeSimulation(t_init);
    gen->V() = 0.0;
    idas->runSimulation(0.1, 20);
    gen->V() = 1.0;
    idas->saveInitialCondition();
  }

  // Get pointer the objective function
  const double* Q = idas->getIntegral();

  // Compute the objective function as an integral over the system trajectory
  idas->getSavedInitialCondition();
  idas->initializeSimulation(t_init);
  idas->initializeQuadrature();
  idas->runSimulationQuadrature(t_final, 100);

  std::cout << "\n\nCost of computing objective function:\n\n";
  idas->printFinalStats();

  const double g1  = Q[0];
  const double eps = 2e-3;

  // Compute gradient of the objective function numerically
  std::vector<double> dGdp(model->sizeParams());

  for (unsigned i = 0; i < model->sizeParams(); ++i)
  {
    model->param()[i] += eps;
    idas->getSavedInitialCondition();
    idas->initializeSimulation(t_init);
    idas->initializeQuadrature();
    idas->runSimulationQuadrature(t_final, 100);

    std::cout << "\n\nCost of computing derivative with respect to parameter "
              << i << ":\n\n";
    idas->printFinalStats();
    double g2 = Q[0];

    // restore parameter to original value
    model->param()[i] -= eps;

    // Evaluate dG/dp numerically
    dGdp[i] = (g2 - g1) / eps;
  }

  // Compute gradient of the objective function using adjoint method
  idas->initializeAdjoint();
  idas->getSavedInitialCondition();
  idas->initializeSimulation(t_init);
  idas->initializeQuadrature();
  idas->runForwardSimulation(t_final, 100);

  std::cout << "\n\nCost of forward simulation for adjoint\n"
            << "sensitivity analysis:\n\n";
  idas->printFinalStats();

  idas->initializeBackwardSimulation(t_final);
  idas->runBackwardSimulation(t_init);

  std::cout << "\n\nCost of adjoint sensitivity analysis:\n\n";
  idas->printFinalStats();

  // Compare results
  int retval = 0;
  std::cout << "\n\nComparison of numerical and adjoint results:\n\n";
  double* neg_dGdp = idas->getAdjointIntegral();
  for (unsigned i = 0; i < model->sizeParams(); ++i)
  {
    std::cout << "dG/dp" << i << " (numerical) = " << dGdp[i] << "\n";
    std::cout << "dG/dp" << i << " (adjoint)   = " << -neg_dGdp[i] << "\n\n";
    if (!isEqual(dGdp[i], -neg_dGdp[i], 10 * eps))
      --retval;
  }

  if (retval < 0)
  {
    std::cout << "The two results differ beyond solver tolerance!\n";
  }

  return retval;
}
