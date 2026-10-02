#include <chrono>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <GridKit/Model/PowerElectronics/SystemModelPowerElectronics.hpp>
#include <GridKit/Solver/Dynamic/DynamicSolver.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>

#include <examples/PowerElectronics/ExamplesHelper/MicrogridNetwork.hpp>
#include <examples/PowerElectronics/ExamplesHelper/SystemAssembler.hpp>

int main(int /* argc */, char const** /* argv */)
{
  double rel_tol         = 1.0e-8;
  size_t max_step_number = 3000;
  bool   use_jac         = true;

  // Create model
  auto* sysmodel = new GridKit::PowerElectronicsModel<double, size_t>(use_jac);

  // Build the four-generator microgrid network.
  size_t                                N_size = 2;
  ScaleMicrogridNetwork<double, size_t> network(N_size);
  assembleSystem(network, *sysmodel);

  // Generator parameters used to construct the initial conditions.
  const auto& parms1 = network.DGParam_list[0];
  const auto& parms2 = network.DGParam_list[2];

  sysmodel->allocate();

  auto* y  = sysmodel->y().getData();
  auto* yp = sysmodel->yp().getData();

  // Create initial points for states.
  for (size_t i = 0; i < sysmodel->size(); ++i)
  {
    y[i]  = 0.0;
    yp[i] = 0.0;
  }

  // Create initial derivatives.
  yp[2]      = parms1.Vn_;
  yp[4]      = parms1.Kpv_ * parms1.Vn_;
  yp[6]      = (parms1.Kpc_ * parms1.Kpv_ * parms1.Vn_) / parms1.Lf_;
  yp[12 + 2] = parms1.Vn_;
  yp[12 + 4] = parms1.Kpv_ * parms1.Vn_;
  yp[12 + 6] = (parms1.Kpc_ * parms1.Kpv_ * parms1.Vn_) / parms1.Lf_;

  for (size_t i = 2; i < 4; ++i)
  {
    yp[13 * i - 1 + 2] = parms2.Vn_;
    yp[13 * i - 1 + 4] = parms2.Kpv_ * parms2.Vn_;
    yp[13 * i - 1 + 6] = (parms2.Kpc_ * parms2.Kpv_ * parms2.Vn_) / parms2.Lf_;
  }

  // Since the initial P_com = 0.
  y[network.dg_signal.getNodeConnection(0).idx_] = parms1.wb_;

  sysmodel->y().setDataUpdated();
  sysmodel->yp().setDataUpdated();

  sysmodel->initialize();
  sysmodel->evaluateResidual();

  sysmodel->updateTime(0.0, 1.0e-8);
  sysmodel->evaluateJacobian();
  sysmodel->tagDifferentiable();

  // ------------------------------------------------------------
  // Output file
  // ------------------------------------------------------------

  std::ofstream output("microgrid_fault.csv");

  output << std::setprecision(16);

  output << "time";
  for (size_t i = 0; i < sysmodel->size(); ++i)
  {
    output << ",y" << i;
  }
  output << '\n';

  auto write_state = [&](double t)
  {
    const auto* state = sysmodel->y().getData();

    output << t;

    for (size_t i = 0; i < sysmodel->size(); ++i)
    {
      output << ',' << state[i];
    }

    output << '\n';
  };

  // ------------------------------------------------------------
  // Create IDA
  // ------------------------------------------------------------

  auto* idas = new AnalysisManager::Sundials::Ida<double, size_t>(sysmodel);

  idas->setTolerance(rel_tol);
  idas->setMaxSteps(max_step_number);
  idas->configureSimulation();
  idas->getDefaultInitialCondition();
  idas->initializeSimulation(0.0);

  // Initial state
  write_state(0.0);

  // ------------------------------------------------------------
  // Run normally: 0 -> 1 s
  // ------------------------------------------------------------

  idas->runSimulation(1.0, 0.001, write_state);

  write_state(1.0);

  std::cout << "Reached t = 1.0 s\n";

  // ------------------------------------------------------------
  // Inject load disturbance at t = 1 s
  // ------------------------------------------------------------

  std::cout << "Injecting load fault at t = 1.0 s\n";

  network.loads[0]->setResistance(3.0);
  network.loads[0]->setInductance(6.4e-3);

  network.loads[2]->setResistance(2.0);
  network.loads[2]->setInductance(3.2e-3);
  const auto t_start = std::chrono::steady_clock::now();

  // ------------------------------------------------------------
  // Reinitialize IDA at the discontinuity.
  //
  // y and yp are NOT reset. They contain the solution reached
  // at t = 1 s.
  // ------------------------------------------------------------
  // idas->initializeSimulation(1.0, true);

  // // ------------------------------------------------------------
  // // Run disturbed system: 1 -> 2 s
  // // ------------------------------------------------------------

  // idas->runSimulation(2.0, 0.01, write_state);

  const auto t_end = std::chrono::steady_clock::now();

  const double elapsed_s = std::chrono::duration<double>(t_end - t_start).count();

  std::cout << "\nCosimulation runtime: " << elapsed_s << " s "
            << " ms/step)\n";

  write_state(2.0);

  std::cout << "Reached t = 2.0 s\n";

  output.close();

  std::cout << "Results written to microgrid_fault.csv\n";

  delete idas;
  delete sysmodel;

  return 0;
}