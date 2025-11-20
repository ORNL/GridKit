

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

int main(int /* argc */, char const** /* argv */)
{
  double abs_tol = 1.0e-8;
  double rel_tol = 1.0e-8;
  bool   use_jac = true;

  // TODO:setup as named parameters
  // Create circuit model
  GridKit::PowerElectronicsModel<double, size_t> sysmodel(rel_tol, abs_tol, use_jac);

  size_t idoff = 0;

  // RL circuit parameters
  double rinit = 1.0;
  double linit = 1.0;
  double vinit = 1.0;

  // System vector map
  // 0 -> Inductor internal
  // 1 -> VSource internal
  // 2 -> Resistor/VSource external
  // 3 -> Inductor/Resistor external

  // inductor
  GridKit::Inductor<double, size_t>* induct = new GridKit::Inductor<double, size_t>(idoff, linit);
  // Form index to node uid realations
  //  input
  induct->setExternalConnectionNodes(0, 3);
  // output
  induct->setExternalConnectionNodes(1, static_cast<size_t>(-1));
  // internal
  induct->setExternalConnectionNodes(2, 0);
  // add component
  sysmodel.addComponent(induct);

  // resistor
  idoff++;
  GridKit::Resistor<double, size_t>* resis = new GridKit::Resistor<double, size_t>(idoff, rinit);
  // Form index to node uid realations
  // input
  resis->setExternalConnectionNodes(0, 2);
  // output
  resis->setExternalConnectionNodes(1, 3);
  // add
  sysmodel.addComponent(resis);

  // voltage source
  idoff++;
  GridKit::VoltageSource<double, size_t>* vsource = new GridKit::VoltageSource<double, size_t>(idoff, vinit);
  // Form index to node uid realations
  // input
  vsource->setExternalConnectionNodes(0, static_cast<size_t>(-1));
  // output
  vsource->setExternalConnectionNodes(1, 2);
  // internal
  vsource->setExternalConnectionNodes(2, 1);

  sysmodel.addComponent(vsource);

  sysmodel.allocate(4);

  std::cout << sysmodel.y().size() << std::endl;

  // Grounding for IDA. If no grounding then circuit is \mu > 1
  // v_0 (grounded)
  // Create Intial points
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

  std::cout << "Verify Intial Resisdual is Zero: {";
  for (double i : sysmodel.getResidual())
  {
    std::cout << i << ", ";
  }
  std::cout << "}\n";

  sysmodel.updateTime(0.0, 1.0);
  sysmodel.evaluateJacobian();
  std::cout << "Intial Jacobian with alpha = 1:\n";
  sysmodel.getJacobian().printMatrix();

  // Create numerical integrator and configure it for the generator model
  AnalysisManager::Sundials::Ida<double, size_t> idas(&sysmodel);

  double t_init  = 0.0;
  double t_final = 1.0;

  // setup simulation
  idas.configureSimulation();
  idas.getDefaultInitialCondition();
  idas.initializeSimulation(t_init);

  idas.runSimulation(t_final);

  std::vector<double>& yfinial = sysmodel.y();

  std::cout << "Final Vector y\n";
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

  std::cout << "Element-wise Relative error at t=" << t_final << "\n";
  for (size_t i = 0; i < yfinial.size(); i++)
  {
    std::cout << abs((yfinial[i] - yexact[i]) / yexact[i]) << "\n";
  }

  return 0;
}
