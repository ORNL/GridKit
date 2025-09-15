

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <GridKit/Model/PowerElectronics/Capacitor/Capacitor.hpp>
#include <GridKit/Model/PowerElectronics/Inductor/Inductor.hpp>
#include <GridKit/Model/PowerElectronics/Resistor/Resistor.hpp>
#include <GridKit/Model/PowerElectronics/SystemModelPowerElectronics.hpp>
#include <GridKit/Model/PowerElectronics/VoltageSource/VoltageSource.hpp>
#include <GridKit/Solver/Dynamic/DynamicSolver.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Utilities/Benchmark.hpp>

struct RLCircuitComponents
{
  std::unique_ptr<GridKit::PowerElectronicsModel<double, size_t>> system_;
  GridKit::Inductor<double, size_t>*                              inductor_;
  GridKit::Resistor<double, size_t>*                              resistor_;
  GridKit::VoltageSource<double, size_t>*                         vsource_;
};

RLCircuitComponents rlCircuit()
{
  double abs_tol = 1.0e-8;
  double rel_tol = 1.0e-8;
  bool   use_jac = true;

  // TODO:setup as named parameters
  // Create circuit model
  auto  sysmodel_ptr = std::make_unique<GridKit::PowerElectronicsModel<double, size_t>>(rel_tol, abs_tol, use_jac);
  auto& sysmodel     = *sysmodel_ptr;

  size_t idoff = 0;

  // RL circuit parameters
  double rinit = 1.0;
  double linit = 1.0;
  double vinit = 1.0;

  // inductor
  auto induct = new GridKit::Inductor<double, size_t>(idoff, linit);
  // Form index to node uid realations
  //  input
  induct->setExternalConnectionNodes(0, 1);
  // output
  induct->setExternalConnectionNodes(1, static_cast<size_t>(-1));
  // internal
  induct->setExternalConnectionNodes(2, 2);
  // add component
  sysmodel.addComponent(induct);

  // resistor
  idoff++;
  auto resis = new GridKit::Resistor<double, size_t>(idoff, rinit);
  // Form index to node uid realations
  // input
  resis->setExternalConnectionNodes(0, 0);
  // output
  resis->setExternalConnectionNodes(1, 1);
  // add
  sysmodel.addComponent(resis);

  // voltage source
  idoff++;
  auto vsource = new GridKit::VoltageSource<double, size_t>(idoff, vinit);
  // Form index to node uid realations
  // input
  vsource->setExternalConnectionNodes(0, static_cast<size_t>(-1));
  // output
  vsource->setExternalConnectionNodes(1, 0);
  // internal
  vsource->setExternalConnectionNodes(2, 3);

  sysmodel.addComponent(vsource);

  sysmodel.allocate(4);

  // Grounding for IDA. If no grounding then circuit is \mu > 1
  // v_0 (grounded)
  // Create Intial points
  sysmodel.y()[0] = vinit; // v_1
  sysmodel.y()[1] = vinit; // v_2
  sysmodel.y()[2] = 0.0;   // i_L
  sysmodel.y()[3] = 0.0;   // i_s

  sysmodel.yp()[0] = 0.0;            // v'_1
  sysmodel.yp()[1] = 0.0;            // v'_2
  sysmodel.yp()[2] = -vinit / linit; // i'_s
  sysmodel.yp()[3] = -vinit / linit; // i'_L

  sysmodel.initialize();

  sysmodel.updateTime(0.0, 1.0);

  // Necessary to be able to calculate the sparsity of the system correctly.
  sysmodel.evaluateJacobian();

  return {.system_ = std::move(sysmodel_ptr), .inductor_ = std::move(induct), .resistor_ = std::move(resis), .vsource_ = std::move(vsource)};
}

void doIterations(bool use_csr, bool csr_inductor, bool csr_resistor, bool csr_vsource)
{
  const size_t NUM_ITERATIONS = 50;

  for (unsigned i = 0; i < NUM_ITERATIONS; i++)
  {
    GridKit::Utility::benchmark.newIteration();

    auto  rl_circuit = rlCircuit();
    auto& sysmodel   = *rl_circuit.system_;

    rl_circuit.inductor_->setCsrUsage(csr_inductor);
    rl_circuit.resistor_->setCsrUsage(csr_resistor);
    rl_circuit.vsource_->setCsrUsage(csr_vsource);

    // Create numerical integrator and configure it for the generator model
    AnalysisManager::Sundials::Ida<double, size_t> idas(&sysmodel);
    idas.setUseCsr(use_csr);

    double t_init  = 0.0;
    double t_final = 1.0;

    // setup simulation
    idas.configureSimulation();
    idas.getDefaultInitialCondition();
    idas.initializeSimulation(t_init);

    auto timer = GridKit::Utility::startTime<true>("RLCircuit Simulation");
    idas.runSimulation(t_final);
    GridKit::Utility::endTime<true>(std::move(timer));
  }
}

int main(int /* argc */, char const** /* argv */)
{
  GridKit::Utility::benchmark.newRun("COO Jacobians");
  doIterations(false, false, false, false);

  GridKit::Utility::benchmark.newRun("CSR - System");
  doIterations(true, false, false, false);

  GridKit::Utility::benchmark.newRun("CSR - Inductor");
  doIterations(true, true, false, false);

  GridKit::Utility::benchmark.newRun("CSR - Inductor, Resistor");
  doIterations(true, true, true, false);

  GridKit::Utility::benchmark.newRun("CSR - All");
  doIterations(true, true, true, true);

  std::cerr << GridKit::Utility::benchmark.report() << '\n';
  return 0;
}