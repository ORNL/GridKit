#pragma once

#include <memory>

#include <Model/PowerElectronics/SystemModelPowerElectronics.hpp>

std::unique_ptr<GridKit::PowerElectronicsModel<double, size_t>> rlCircuitSystem(double abs_tol, double rel_tol, bool use_jac, double linit, double rinit, double vinit)
{
  // TODO:setup as named parameters
  // Create circuit model
  auto  sysmodel_ptr = std::make_unique<GridKit::PowerElectronicsModel<double, size_t>>(rel_tol, abs_tol, use_jac);
  auto& sysmodel     = *sysmodel_ptr;

  size_t idoff = 0;

  // inductor
  GridKit::Inductor<double, size_t>* induct = new GridKit::Inductor<double, size_t>(idoff, linit);
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
  GridKit::Resistor<double, size_t>* resis = new GridKit::Resistor<double, size_t>(idoff, rinit);
  // Form index to node uid realations
  // input
  resis->setExternalConnectionNodes(0, 0);
  // output
  resis->setExternalConnectionNodes(1, 1);
  // add
  sysmodel.addComponent(resis);

  // voltage source
  idoff++;
  GridKit::VoltageSource<double, size_t>* vsource = new GridKit::VoltageSource<double, size_t>(idoff, vinit);
  // Form index to node uid realations
  // input
  vsource->setExternalConnectionNodes(0, static_cast<size_t>(-1));
  // output
  vsource->setExternalConnectionNodes(1, 0);
  // internal
  vsource->setExternalConnectionNodes(2, 3);

  sysmodel.addComponent(vsource);

  sysmodel.allocate(4);

  std::cout << sysmodel.y().size() << std::endl;

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

  return sysmodel_ptr;
}
