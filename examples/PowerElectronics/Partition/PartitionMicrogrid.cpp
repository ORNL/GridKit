#include <algorithm>
#include <cstddef>

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>

#include <GridKit/Model/PowerElectronics/DistributedGenerator/DistributedGenerator.hpp>
#include <GridKit/Model/PowerElectronics/MicrogridBusDQ/MicrogridBusDQ.hpp>
#include <GridKit/Model/PowerElectronics/MicrogridLine/MicrogridLine.hpp>
#include <GridKit/Model/PowerElectronics/MicrogridLoad/MicrogridLoad.hpp>
#include <GridKit/Model/PowerElectronics/PartitionInterface/BusPartitionInterface.hpp>
#include <GridKit/Model/PowerElectronics/SubsystemModel.hpp>
#include <GridKit/Model/PowerElectronics/SystemModelPowerElectronics.hpp>

int main()
{
  /// @todo Needs to be modified. Some components are small relative to others thus
  /// there error is high (or could be matlab vector issue)
  double abs_tol         = 1.0e-8;
  double rel_tol         = 1.0e-8;
  size_t max_step_number = 3000;
  bool   use_jac         = true;

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

  // indexing sets
  size_t Nsize              = 2;
  //							DGs	+		- refframe	   Lines +				Loads
  size_t vec_size_internals = 13 * (2 * Nsize) - 1 + (2 + 4 * (Nsize - 1)) + 2 * Nsize;
  //							\omegaref + BusDQ
  size_t vec_size_externals = 1 + 2 * (2 * Nsize);
  size_t dqbus1             = vec_size_internals + 1;
  size_t dqbus2             = vec_size_internals + 3;
  size_t dqbus3             = vec_size_internals + 5;
  size_t dqbus4             = vec_size_internals + 7;

  size_t vec_size_total = vec_size_internals + vec_size_externals;

  size_t indexv = 0;

  // dg 1
  GridKit::DistributedGenerator<double, size_t>* dg1 = new GridKit::DistributedGenerator<double, size_t>(0, parms1, true);
  // ref motor
  dg1->setExternalConnectionNodes(0, vec_size_internals);
  // outputs
  dg1->setExternalConnectionNodes(1, dqbus1);
  dg1->setExternalConnectionNodes(2, dqbus1 + 1);
  //"grounding" of the difference
  dg1->setExternalConnectionNodes(3, static_cast<size_t>(-1));
  // internal connections
  for (size_t i = 0; i < 12; i++)
  {

    dg1->setExternalConnectionNodes(4 + i, indexv + i);
  }
  indexv += 12;

  // dg 2
  GridKit::DistributedGenerator<double, size_t>* dg2 = new GridKit::DistributedGenerator<double, size_t>(1, parms1, false);
  // ref motor
  dg2->setExternalConnectionNodes(0, vec_size_internals);
  // outputs
  dg2->setExternalConnectionNodes(1, dqbus2);
  dg2->setExternalConnectionNodes(2, dqbus2 + 1);
  // internal connections
  for (size_t i = 0; i < 13; i++)
  {

    dg2->setExternalConnectionNodes(3 + i, indexv + i);
  }
  indexv += 13;

  // dg 3
  GridKit::DistributedGenerator<double, size_t>* dg3 = new GridKit::DistributedGenerator<double, size_t>(2, parms2, false);
  // ref motor
  dg3->setExternalConnectionNodes(0, vec_size_internals);
  // outputs
  dg3->setExternalConnectionNodes(1, dqbus3);
  dg3->setExternalConnectionNodes(2, dqbus3 + 1);
  // internal connections
  for (size_t i = 0; i < 13; i++)
  {

    dg3->setExternalConnectionNodes(3 + i, indexv + i);
  }
  indexv += 13;

  // dg 4
  GridKit::DistributedGenerator<double, size_t>* dg4 = new GridKit::DistributedGenerator<double, size_t>(3, parms2, false);
  // ref motor
  dg4->setExternalConnectionNodes(0, vec_size_internals);
  // outputs
  dg4->setExternalConnectionNodes(1, dqbus4);
  dg4->setExternalConnectionNodes(2, dqbus4 + 1);

  // internal connections
  for (size_t i = 0; i < 13; i++)
  {

    dg4->setExternalConnectionNodes(3 + i, indexv + i);
  }
  indexv += 13;

  // Lines

  // line 1
  GridKit::MicrogridLine<double, size_t>* l1 = new GridKit::MicrogridLine<double, size_t>(4, rline1, Lline1);
  // ref motor
  l1->setExternalConnectionNodes(0, vec_size_internals);
  // input connections
  l1->setExternalConnectionNodes(1, dqbus1);
  l1->setExternalConnectionNodes(2, dqbus1 + 1);
  // output connections
  l1->setExternalConnectionNodes(3, dqbus2);
  l1->setExternalConnectionNodes(4, dqbus2 + 1);
  // internal connections
  for (size_t i = 0; i < 2; i++)
  {

    l1->setExternalConnectionNodes(5 + i, indexv + i);
  }
  indexv += 2;

  // line 2
  GridKit::MicrogridLine<double, size_t>* l2 = new GridKit::MicrogridLine<double, size_t>(5, rline2, Lline2);
  // ref motor
  l2->setExternalConnectionNodes(0, vec_size_internals);
  // input connections
  l2->setExternalConnectionNodes(1, dqbus2);
  l2->setExternalConnectionNodes(2, dqbus2 + 1);
  // output connections
  l2->setExternalConnectionNodes(3, dqbus3);
  l2->setExternalConnectionNodes(4, dqbus3 + 1);
  // internal connections
  for (size_t i = 0; i < 2; i++)
  {

    l2->setExternalConnectionNodes(5 + i, indexv + i);
  }
  indexv += 2;

  // line 3
  GridKit::MicrogridLine<double, size_t>* l3 = new GridKit::MicrogridLine<double, size_t>(6, rline3, Lline3);
  // ref motor
  l3->setExternalConnectionNodes(0, vec_size_internals);
  // input connections
  l3->setExternalConnectionNodes(1, dqbus3);
  l3->setExternalConnectionNodes(2, dqbus3 + 1);
  // output connections
  l3->setExternalConnectionNodes(3, dqbus4);
  l3->setExternalConnectionNodes(4, dqbus4 + 1);
  // internal connections
  for (size_t i = 0; i < 2; i++)
  {

    l3->setExternalConnectionNodes(5 + i, indexv + i);
  }
  indexv += 2;

  //  loads

  // load 1
  GridKit::MicrogridLoad<double, size_t>* load1 = new GridKit::MicrogridLoad<double, size_t>(7, rload1, Lload1);
  // ref motor
  load1->setExternalConnectionNodes(0, vec_size_internals);
  // input connections
  load1->setExternalConnectionNodes(1, dqbus1);
  load1->setExternalConnectionNodes(2, dqbus1 + 1);
  // internal connections
  for (size_t i = 0; i < 2; i++)
  {

    load1->setExternalConnectionNodes(3 + i, indexv + i);
  }
  indexv += 2;

  // load 2
  GridKit::MicrogridLoad<double, size_t>* load2 = new GridKit::MicrogridLoad<double, size_t>(8, rload2, Lload2);
  // ref motor
  load2->setExternalConnectionNodes(0, vec_size_internals);
  // input connections
  load2->setExternalConnectionNodes(1, dqbus3);
  load2->setExternalConnectionNodes(2, dqbus3 + 1);
  // internal connections
  for (size_t i = 0; i < 2; i++)
  {

    load2->setExternalConnectionNodes(3 + i, indexv + i);
  }
  indexv += 2;

  // Virtual PQ Buses
  GridKit::MicrogridBusDQ<double, size_t>* bus1 = new GridKit::MicrogridBusDQ<double, size_t>(9, RN);

  bus1->setExternalConnectionNodes(0, dqbus1);
  bus1->setExternalConnectionNodes(1, dqbus1 + 1);

  GridKit::MicrogridBusDQ<double, size_t>* bus2 = new GridKit::MicrogridBusDQ<double, size_t>(10, RN);

  bus2->setExternalConnectionNodes(0, dqbus2);
  bus2->setExternalConnectionNodes(1, dqbus2 + 1);

  GridKit::MicrogridBusDQ<double, size_t>* bus3 = new GridKit::MicrogridBusDQ<double, size_t>(11, RN);

  bus3->setExternalConnectionNodes(0, dqbus3);
  bus3->setExternalConnectionNodes(1, dqbus3 + 1);

  GridKit::MicrogridBusDQ<double, size_t>* bus4 = new GridKit::MicrogridBusDQ<double, size_t>(12, RN);

  bus4->setExternalConnectionNodes(0, dqbus4);
  bus4->setExternalConnectionNodes(1, dqbus4 + 1);

  GridKit::SubsystemModel<double, size_t>* partition1 = new GridKit::SubsystemModel<double, size_t>(vec_size_internals);
  GridKit::SubsystemModel<double, size_t>* partition2 = new GridKit::SubsystemModel<double, size_t>();

  GridKit::MicrogridLine<double, size_t>          comp3copy(*l2);
  GridKit::BusPartitionInterface<double, size_t>* busInterface = new GridKit::BusPartitionInterface<double, size_t>(comp3copy, dqbus2, dqbus2 + 1, 14);

  partition1->addComponent(dg1);
  partition1->addComponent(dg2);
  partition1->addComponent(l1);
  partition1->addComponent(bus1);
  partition1->addComponent(bus2);
  partition1->addComponent(load1);
  partition1->addComponent(busInterface);

  partition2->addComponent(dg3);
  partition2->addComponent(dg4);
  partition2->addComponent(l2);
  partition2->addComponent(bus3);
  partition2->addComponent(bus4);
  partition2->addComponent(l3);
  partition2->addComponent(load2);

  sysmodel->addComponent(new GridKit::DistributedGenerator<double, size_t>(*dg1));
  sysmodel->addComponent(new GridKit::DistributedGenerator<double, size_t>(*dg2));
  sysmodel->addComponent(new GridKit::DistributedGenerator<double, size_t>(*dg3));
  sysmodel->addComponent(new GridKit::DistributedGenerator<double, size_t>(*dg4));
  sysmodel->addComponent(new GridKit::MicrogridLine<double, size_t>(*l1));
  sysmodel->addComponent(new GridKit::MicrogridLine<double, size_t>(*l2));
  sysmodel->addComponent(new GridKit::MicrogridLine<double, size_t>(*l3));
  sysmodel->addComponent(new GridKit::MicrogridBusDQ<double, size_t>(*bus1));
  sysmodel->addComponent(new GridKit::MicrogridBusDQ<double, size_t>(*bus2));
  sysmodel->addComponent(new GridKit::MicrogridBusDQ<double, size_t>(*bus3));
  sysmodel->addComponent(new GridKit::MicrogridBusDQ<double, size_t>(*bus4));
  sysmodel->addComponent(new GridKit::MicrogridLoad<double, size_t>(*load1));
  sysmodel->addComponent(new GridKit::MicrogridLoad<double, size_t>(*load2));

  partition1->allocate();
  partition2->allocate();
  sysmodel->allocate(vec_size_total);

  std::vector<double> y;
  std::vector<double> yp;

  for (size_t i = 0; i < vec_size_total; i++)
  {
    y.push_back(static_cast<double>(i + 1));
    yp.push_back(static_cast<double>(i + 1));
  }

  // Distribute externals to partition 1
  for (size_t i = 0; i < partition1->getExternSize(); i++)
  {
    partition1->getExternalDataY()[i]  = y[partition1->getExternalIndices()[i]];
    partition1->getExternalDataYP()[i] = yp[partition1->getExternalIndices()[i]];
  }

  // Distribute externals to partition 2
  for (size_t i = 0; i < partition2->getExternSize(); i++)
  {
    partition2->getExternalDataY()[i]  = y[partition2->getExternalIndices()[i]];
    partition2->getExternalDataYP()[i] = yp[partition2->getExternalIndices()[i]];
  }

  // Distribute internals to partition 1
  for (size_t i = 0; i < partition1->getInternalSize(); i++)
  {
    partition1->y()[i]  = y[partition1->getNodeConnection(i)];
    partition1->yp()[i] = yp[partition1->getNodeConnection(i)];
  }

  // Distribute internals to partition 2
  for (size_t i = 0; i < partition2->getInternalSize(); i++)
  {
    partition2->y()[i]  = y[partition2->getNodeConnection(i)];
    partition2->yp()[i] = yp[partition2->getNodeConnection(i)];
  }

  // Evaluate Residuals for each partition
  partition1->evaluateResidual();
  partition2->evaluateResidual();

  auto printTitle = [](std::string msg) -> void
  {
    std::cout << "\n--------------- " << msg << " -------------" << std::endl;
    printf("%-12s  ----------  %12s\n", "Res Values", "Comp. Index");
  };

  // Print Residuals from partition 1
  printTitle("Partition 1");
  for (size_t i = 0; i < partition1->getInternalSize(); i++)
  {
    auto com_index = partition1->getNodeConnection(static_cast<size_t>(i));
    printf("%-12.5g  ----------  %7zu\n", partition1->getResidual()[i], com_index);
  }

  // Print Residuals from partition 2
  printTitle("Partition 2");
  for (size_t i = 0; i < partition2->getInternalSize(); i++)
  {
    auto com_index = partition2->getNodeConnection(static_cast<size_t>(i));
    printf("%-12.5g  ----------  %7zu\n", partition2->getResidual()[i], com_index);
  }

  for (size_t i = 0; i < vec_size_total; i++)
  {
    sysmodel->y()[i]  = y[i];
    sysmodel->yp()[i] = yp[i];
  }

  sysmodel->evaluateResidual();

  printTitle("Reference Solution");
  for (size_t i = 0; i < vec_size_total; i++)
  {
    printf("%-12.5g  ----------  %7zu\n", sysmodel->getResidual()[i], i);
  }

  std::vector<double> f(vec_size_total, 0.0);
  std::vector<double> error(vec_size_total, 1.0);

  // Get internal residuals from partition 1
  for (size_t i = 0; i < partition1->getInternalSize(); i++)
  {
    f[partition1->getNodeConnection(i)] = partition1->getResidual()[i];
  }

  // Get internal residuals from partition 2
  for (size_t i = 0; i < partition2->getInternalSize(); i++)
  {
    f[partition2->getNodeConnection(i)] = partition2->getResidual()[i];
  }

  for (size_t i = 0; i < vec_size_total; i++)
  {
    error[i] = sysmodel->getResidual()[i] - f[i];
  }

  double max_error = 0;
  for (size_t i = 0; i < vec_size_total; i++)
  {
    if (max_error < std::abs(error[i]))
    {
      max_error = std::abs(error[i]);
    }
  }

  std::cout << "\nMax Error of Reference and Partition Evaluation: " << max_error << std::endl;

  delete partition1;
  delete partition2;
  delete sysmodel;

  return 0;
}