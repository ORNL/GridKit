

#include "RLCircuit.hpp"

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

  // RL circuit parameters
  double rinit = 1.0;
  double linit = 1.0;
  double vinit = 1.0;

  auto  sysmodel_ptr = rlCircuitSystem(abs_tol, rel_tol, use_jac, linit, rinit, vinit);
  auto& sysmodel     = *sysmodel_ptr;

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
  yexact[0] = vinit;
  yexact[2] = (vinit / rinit) * (exp(-(rinit / linit) * t_final) - 1.0);
  yexact[3] = yexact[2];
  yexact[1] = vinit + rinit * yexact[2];

  std::cout << "Element-wise Relative error at t=" << t_final << "\n";
  for (size_t i = 0; i < yfinial.size(); i++)
  {
    std::cout << abs((yfinial[i] - yexact[i]) / yexact[i]) << "\n";
  }

  return 0;
}
