

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>

#include <GridKit/Model/PowerElectronics/Bus/Bus.hpp>
#include <GridKit/Model/PowerElectronics/Bus/GroundedBus.hpp>
#include <GridKit/Model/PowerElectronics/Capacitor/Capacitor.hpp>
#include <GridKit/Model/PowerElectronics/Inductor/Inductor.hpp>
#include <GridKit/Model/PowerElectronics/Resistor/Resistor.hpp>
#include <GridKit/Model/PowerElectronics/SystemModelPowerElectronics.hpp>
#include <GridKit/Model/PowerElectronics/VoltageSource/VoltageSource.hpp>
#include <GridKit/Solver/Dynamic/DynamicSolver.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>

int main(int /* argc */, char const** /* argv */)
{
  double abs_tol = 1.0e-8;
  double rel_tol = 1.0e-8;
  bool   use_jac = true;

  // TODO:setup as named parameters
  // Create circuit model
  GridKit::PowerElectronicsModel<double, size_t> sysmodel(use_jac);

  size_t idoff = 0;

  // RL circuit parameters
  double rinit = 1.0;
  double linit = 1.0;
  double vinit = 1.0;

  using Bus         = GridKit::PowerElectronics::Bus<double, size_t>;
  using GroundedBus = GridKit::PowerElectronics::GroundedBus<double, size_t>;
  GroundedBus bus_iv(0.0);
  Bus         bus_vr;
  Bus         bus_ir;

  sysmodel.addNode(&bus_iv);
  sysmodel.addNode(&bus_vr);
  sysmodel.addNode(&bus_ir);

  // inductor
  GridKit::Inductor<double, size_t>* induct = new GridKit::Inductor<double, size_t>(idoff, linit, &bus_ir, &bus_iv);
  sysmodel.addComponent(induct);

  // resistor
  idoff++;
  GridKit::Resistor<double, size_t>* resis = new GridKit::Resistor<double, size_t>(idoff, rinit, &bus_vr, &bus_ir);
  sysmodel.addComponent(resis);

  // voltage source
  idoff++;
  GridKit::VoltageSource<double, size_t>* vsource = new GridKit::VoltageSource<double, size_t>(idoff, vinit, &bus_iv, &bus_vr);
  sysmodel.addComponent(vsource);

  sysmodel.allocate();

  std::cout << sysmodel.y().size() << std::endl;

  // Grounding for IDA. If no grounding then circuit is \mu > 1
  // v_0 (grounded)
  // Create initial points
  sysmodel.y()[0] = 0.0;   // i_L
  sysmodel.y()[1] = 0.0;   // i_s
  sysmodel.y()[2] = vinit; // v_1
  sysmodel.y()[3] = vinit; // v_2

  sysmodel.yp()[0] = -vinit / linit; // i'_s
  sysmodel.yp()[1] = -vinit / linit; // i'_L
  sysmodel.yp()[2] = 0.0;            // v'_1
  sysmodel.yp()[3] = 0.0;            // v'_2

  sysmodel.initialize();
  sysmodel.evaluateResidual();

  std::cout << "Verify initial resisdual is zero: {";
  for (double i : sysmodel.getResidual())
  {
    std::cout << i << ", ";
  }
  std::cout << "}\n";

  sysmodel.updateTime(0.0, 1.0);
  sysmodel.evaluateJacobian();
  std::cout << "Initial Jacobian with alpha = 1:\n";
  sysmodel.getCsrJacobian()->print();

  // Create numerical integrator and configure it for the generator model
  AnalysisManager::Sundials::Ida<double, size_t> idas(&sysmodel);

  double t_init  = 0.0;
  double t_final = 1.0;

  // setup simulation
  idas.configureSimulation();
  idas.setTolerance(rel_tol, abs_tol);
  idas.getDefaultInitialCondition();
  idas.initializeSimulation(t_init);

  idas.runSimulation(t_final);

  std::vector<double>& yfinial = sysmodel.y();

  std::cout << "Final vector y\n";
  for (size_t i = 0; i < yfinial.size(); i++)
  {
    std::cout << yfinial[i] << "\n";
  }

  std::vector<double> yexact(4);

  // analytical solution to the circuit
  yexact[2] = vinit;
  yexact[0] = (vinit / rinit) * (exp(-(rinit / linit) * t_final) - 1.0);
  yexact[1] = yexact[0];
  yexact[3] = vinit + rinit * yexact[0];

  std::cout << "Element-wise relative error at t=" << t_final << "\n";
  for (size_t i = 0; i < yfinial.size(); i++)
  {
    std::cout << abs((yfinial[i] - yexact[i]) / yexact[i]) << "\n";
  }

  std::cerr << idas.getStats().report() << '\n';

  return 0;
}
