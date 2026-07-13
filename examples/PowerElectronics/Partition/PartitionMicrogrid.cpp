
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>

#include <GridKit/Model/PowerElectronics/Bus/MicrogridBus.hpp>
#include <GridKit/Model/PowerElectronics/Bus/SignalNode.hpp>
#include <GridKit/Model/PowerElectronics/DistributedGenerator/DistributedGenerator.hpp>
#include <GridKit/Model/PowerElectronics/MicrogridBusDQ/MicrogridBusDQ.hpp>
#include <GridKit/Model/PowerElectronics/MicrogridLine/MicrogridLine.hpp>
#include <GridKit/Model/PowerElectronics/MicrogridLoad/MicrogridLoad.hpp>
#include <GridKit/Model/PowerElectronics/PartitionInterface/BusPartitionInterface.hpp>
#include <GridKit/Model/PowerElectronics/SubsystemModel.hpp>
#include <GridKit/Model/PowerElectronics/SystemModelPowerElectronics.hpp>
#include <GridKit/Solver/Dynamic/DynamicSolver.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>

#include "jac_test_helper.hpp"

int main()
{
  /// @todo Needs to be modified. Some components are small relative to others thus
  /// there error is high (or could be matlab vector issue)
  double abs_tol         = 1.0e-8;
  double rel_tol         = 1.0e-8;
  size_t max_step_number = 3000;
  bool   use_jac         = true;

  // Create model
  auto* sysmodel = new GridKit::PowerElectronicsModel<double, size_t>(rel_tol, abs_tol, use_jac, max_step_number);

  // Modeled after the problem in the paper
  double RN = 1.0e4;

  // DG Params
  GridKit::DistributedGeneratorParameters<double, size_t> parms1;
  parms1.wb_  = 2.0 * M_PI * 50.0;
  parms1.wc_  = 31.41;
  parms1.mp_  = 9.4e-5;
  parms1.Vn_  = 380.0;
  parms1.nq_  = 1.3e-3;
  parms1.F_   = 0.75;
  parms1.Kiv_ = 420.0;
  parms1.Kpv_ = 0.1;
  parms1.Kic_ = 2.0e4;
  parms1.Kpc_ = 15.0;
  parms1.Cf_  = 5.0e-5;
  parms1.rLf_ = 0.1;
  parms1.Lf_  = 1.35e-3;
  parms1.rLc_ = 0.03;
  parms1.Lc_  = 0.35e-3;

  GridKit::DistributedGeneratorParameters<double, size_t> parms2;
  // Parameters from MATLAB Microgrid code for first DG
  parms2.wb_  = 2.0 * M_PI * 50.0;
  parms2.wc_  = 31.41;
  parms2.mp_  = 12.5e-5;
  parms2.Vn_  = 380.0;
  parms2.nq_  = 1.5e-3;
  parms2.F_   = 0.75;
  parms2.Kiv_ = 390.0;
  parms2.Kpv_ = 0.05;
  parms2.Kic_ = 16.0e3;
  parms2.Kpc_ = 10.5;
  parms2.Cf_  = 50.0e-6;
  parms2.rLf_ = 0.1;
  parms2.Lf_  = 1.35e-3;
  parms2.rLc_ = 0.03;
  parms2.Lc_  = 0.35e-3;

  // Line params
  double rline1 = 0.23;
  double Lline1 = 0.1 / (2.0 * M_PI * 50.0);

  double rline2 = 0.35;
  double Lline2 = 0.58 / (2.0 * M_PI * 50.0);

  double rline3 = 0.23;
  double Lline3 = 0.1 / (2.0 * M_PI * 50.0);

  // load parms
  double rload1 = 3.0;
  double Lload1 = 2.0 / (2.0 * M_PI * 50.0);

  double rload2 = 2.0;
  double Lload2 = 1.0 / (2.0 * M_PI * 50.0);

  using SignalNode = GridKit::PowerElectronics::SignalNode<double, size_t>;
  SignalNode dg_signal;

  sysmodel->addNode(&dg_signal);

  using Bus = GridKit::PowerElectronics::MicrogridBus<double, size_t>;
  Bus bus1;
  Bus bus2;
  Bus bus3;
  Bus bus4;

  sysmodel->addNode(&bus1);
  sysmodel->addNode(&bus2);
  sysmodel->addNode(&bus3);
  sysmodel->addNode(&bus4);

  // dg 1
  GridKit::DistributedGenerator<double, size_t>* dg1 = new GridKit::DistributedGenerator<double, size_t>(
      0, parms1, true, &dg_signal, &bus1);
  sysmodel->addComponent(dg1);

  // dg 2
  GridKit::DistributedGenerator<double, size_t>* dg2 = new GridKit::DistributedGenerator<double, size_t>(
      1, parms1, false, &dg_signal, &bus2);
  sysmodel->addComponent(dg2);

  // dg 3
  GridKit::DistributedGenerator<double, size_t>* dg3 = new GridKit::DistributedGenerator<double, size_t>(
      2, parms2, false, &dg_signal, &bus3);
  sysmodel->addComponent(dg3);

  // dg 4
  GridKit::DistributedGenerator<double, size_t>* dg4 = new GridKit::DistributedGenerator<double, size_t>(
      3, parms2, false, &dg_signal, &bus4);
  sysmodel->addComponent(dg4);

  // Lines

  // line 1
  GridKit::MicrogridLine<double, size_t>* l1 = new GridKit::MicrogridLine<double, size_t>(
      4, rline1, Lline1, &dg_signal, &bus1, &bus2);
  sysmodel->addComponent(l1);

  // line 2
  GridKit::MicrogridLine<double, size_t>* l2 = new GridKit::MicrogridLine<double, size_t>(
      5, rline2, Lline2, &dg_signal, &bus2, &bus3);
  sysmodel->addComponent(l2);

  // line 3
  GridKit::MicrogridLine<double, size_t>* l3 = new GridKit::MicrogridLine<double, size_t>(
      6, rline3, Lline3, &dg_signal, &bus3, &bus4);
  sysmodel->addComponent(l3);

  //  loads

  // load 1
  GridKit::MicrogridLoad<double, size_t>* load1 = new GridKit::MicrogridLoad<double, size_t>(7, rload1, Lload1, &dg_signal, &bus1);
  sysmodel->addComponent(load1);

  // load 2
  GridKit::MicrogridLoad<double, size_t>* load2 = new GridKit::MicrogridLoad<double, size_t>(8, rload2, Lload2, &dg_signal, &bus3);
  sysmodel->addComponent(load2);

  // Virtual PQ Buses
  GridKit::MicrogridBusDQ<double, size_t>* bus_para_1 = new GridKit::MicrogridBusDQ<double, size_t>(9, RN, &bus1);
  sysmodel->addComponent(bus_para_1);

  GridKit::MicrogridBusDQ<double, size_t>* bus_para_2 = new GridKit::MicrogridBusDQ<double, size_t>(10, RN, &bus2);
  sysmodel->addComponent(bus_para_2);

  GridKit::MicrogridBusDQ<double, size_t>* bus_para_3 = new GridKit::MicrogridBusDQ<double, size_t>(11, RN, &bus3);
  sysmodel->addComponent(bus_para_3);

  GridKit::MicrogridBusDQ<double, size_t>* bus_para_4 = new GridKit::MicrogridBusDQ<double, size_t>(12, RN, &bus4);
  sysmodel->addComponent(bus_para_4);

  sysmodel->allocate();

  GridKit::SubsystemModel<double, size_t>* partition1 = new GridKit::SubsystemModel<double, size_t>();
  GridKit::SubsystemModel<double, size_t>* partition2 = new GridKit::SubsystemModel<double, size_t>();

  GridKit::MicrogridLine<double, size_t>          l2copy(*l2);
  GridKit::BusPartitionInterface<double, size_t>* busInterface1 = new GridKit::BusPartitionInterface<double, size_t>(bus2, l2copy, 14);

  busInterface1->allocate();

  partition1->addNode(&dg_signal);
  partition1->addComponent(dg1);
  partition1->addComponent(dg2);
  partition1->addComponent(l1);
  partition1->addComponent(load1);
  partition1->addComponent(busInterface1);
  partition1->addComponent(bus_para_1);
  partition1->addComponent(bus_para_2);
  partition1->addNode(&bus1);
  partition1->addNode(&bus2);

  partition2->addComponent(dg3);
  partition2->addComponent(dg4);
  partition2->addComponent(l2);
  partition2->addComponent(l3);
  partition2->addComponent(load2);
  partition2->addComponent(bus_para_4);
  partition2->addComponent(bus_para_3);
  partition2->addNode(&bus3);
  partition2->addNode(&bus4);

  std::vector<double> y;
  std::vector<double> yp;

  for (size_t i = 0; i < sysmodel->size(); i++)
  {
    y.push_back(static_cast<double>(i + 1));
    yp.push_back(static_cast<double>(i + 1));
  }

  for (size_t i = 0; i < sysmodel->size(); i++)
  {
    sysmodel->y()[i]  = y[i];
    sysmodel->yp()[i] = yp[i];
  }

  sysmodel->updateTime(2, 5);
  sysmodel->evaluateResidual();
  sysmodel->evaluateJacobian();

  auto full_jac = sysmodel->getCsrJacobian();

  std::vector<double> f_sysmodel = sysmodel->getResidual();

  std::vector<GridKit::SubsystemModel<double, size_t>*> partitions = {partition1, partition2};

  // Distribute externals to partition 1
  for (auto* partition : partitions)
  {
    partition->allocate();

    partition->updateTime(2, 5);

    for (size_t i = 0; i < partition->getExternSize(); i++)
    {
      partition->getExternalDataY()[i]  = y[partition->getExternalIndices()[i]];
      partition->getExternalDataYP()[i] = yp[partition->getExternalIndices()[i]];
    }

    for (size_t i = 0; i < partition->getInternalSize(); i++)
    {
      partition->y()[i]  = y[partition->getNodeConnection(i)];
      partition->yp()[i] = yp[partition->getNodeConnection(i)];
    }

    partition->evaluateResidual();
    partition->evaluateJacobian();
  }

  auto partition1_jac = partition1->getCsrJacobian();
  auto partition2_jac = partition2->getCsrJacobian();

  bool jac_matched_1 = GridKit::Testing::verifySubsystemJacobian(*full_jac, *partition1_jac, *partition1);
  bool jac_matched_2 = GridKit::Testing::verifySubsystemJacobian(*full_jac, *partition2_jac, *partition2);

  std::cout << "Jacobian Matched: " << jac_matched_1 * jac_matched_2 << std::endl;

  if (!jac_matched_1 || !jac_matched_2)
  {
    std::cout << "ERROR: At least one subsystem Jacobian is incorrect!" << std::endl;
    return 1;
  }

  std::vector<double> f(sysmodel->size(), 0.0);
  std::vector<double> error(sysmodel->size(), 1.0);

  // Get internal residuals from partition 1
  for (auto* partition : partitions)
  {
    for (size_t i = 0; i < partition->getInternalSize(); i++)
    {
      f[partition->getNodeConnection(i)] = partition->getResidual()[i];
    }
  }

  double max_error = 0;
  for (size_t i = 0; i < sysmodel->size(); i++)
  {
    error[i] = (f_sysmodel[i] - f[i]) / f_sysmodel[i];
    if (max_error < std::abs(error[i]))
    {
      max_error = std::abs(error[i]);
    }
  }

  if (max_error > std::numeric_limits<double>::epsilon())
  {
    std::cout << "ERROR: Max Error too high!" << std::endl;
    return 1;
  }

  std::cout << "\nMax Error of Reference and Partition Evaluation: " << max_error << std::endl;

  delete sysmodel;
  delete partition1;
  delete partition2;

  return 0;
}
