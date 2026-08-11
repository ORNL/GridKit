#include <cmath>
#include <cstddef>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <string>

#include <GridKit/Model/PowerElectronics/Bus/MicrogridBus.hpp>
#include <GridKit/Model/PowerElectronics/PartitionInterface/BusPartitionInterface.hpp>
#include <GridKit/Model/PowerElectronics/SubsystemModel.hpp>
#include <GridKit/Model/PowerElectronics/SystemModelPowerElectronics.hpp>

#include "HiresBus.hpp"
#include "HiresComponent1.hpp"
#include "HiresComponent3.hpp"

std::vector<double> hires(const std::vector<double>& y, const std::vector<double>& yp);
void                printRow(double residual, size_t index);
void                printHeader(const std::string& title);
void                printPartitionResiduals(auto* partition, const std::string& title);

/**
 * This example assembles the HIRES problem and partitions it to demonstrate
 * that the partitioned system can reconstruct the full system.
 *
 * The HIRES problem consists of eight ODE equations and is divided into three
 * components. HiresComponent1 contains equations 1, 2, and 3,
 * HiresComponent3 contains equations 6, 7, and 8, and HiresBus contains
 * equations 4 and 5. HiresBus plays a role similar to MicrogridBusDQ, where
 * other components connected to the bus contribute terms to its equations.
 *
 * In the equations below, terms enclosed in parentheses ( ) are contributions
 * from HiresComponent1, terms enclosed in square brackets [ ] are contributions
 * from HiresComponent3, and the remaining terms in equations 4 and 5 belong
 * to HiresBus.
 *
 * f_1 = dy_1/dt + 1.71y_1 - 0.43y_2 - 8.32y_3 - 0.0007
 *
 * f_2 = dy_2/dt - 1.71y_1 + 8.75y_2
 *
 * f_3 = dy_3/dt + 10.03y_3 - 0.43y_4 - 0.035y_5
 *
 * f_4 = dy_4/dt + y_4 + (0.1y_4 - 8.32y_2 - 1.71y_3) + [0.02y_4]
 *
 * f_5 = dy_5/dt + y_5 + (0.7y_5) + [0.045y_5 - 0.43y_6 - 0.43y_7]
 *
 * f_6 = dy_6/dt - 280y_6y_8 + 0.69y_4 + 1.71y_5 - 0.43y_6 + 0.69y_7
 *
 * f_7 = dy_7/dt + 280y_6y_8 - 1.81y_7
 *
 * f_8 = dy_8/dt - 280y_6y_8 + 1.81y_7
 *
 *
 * The assembled system has the following structure:
 *
 * (HiresComponent1) -------- (HiresBus) -------- (HiresComponent3)
 *
 * The system is partitioned between HiresBus and HiresComponent3. A bus
 * partition interface is introduced in the partition containing HiresBus to
 * preserve the contribution of HiresComponent3 to the bus equations. The
 * residuals evaluated independently by the partitions are then printed with
 * the residual of the full system for eye-ball comparision.
 *
 * This example will also be useful later for testing the order of accuracy of co-simulation methods.
 */

int main(int /* argc */, char const** /* argv */)
{

  using Bus = GridKit::PowerElectronics::MicrogridBus<double, size_t>;
  Bus bus;

  GridKit::PowerElectronicsModel<double, size_t>* sys        = new GridKit::PowerElectronicsModel<double, size_t>();
  GridKit::SubsystemModel<double, size_t>*        partition1 = new GridKit::SubsystemModel<double, size_t>();
  GridKit::SubsystemModel<double, size_t>*        partition2 = new GridKit::SubsystemModel<double, size_t>();

  GridKit::HiresComponent1<double, size_t>* comp1 = new GridKit::HiresComponent1<double, size_t>(&bus, 1);
  GridKit::HiresBus<double, size_t>*        bus1  = new GridKit::HiresBus<double, size_t>(&bus, 2);
  GridKit::HiresComponent3<double, size_t>* comp3 = new GridKit::HiresComponent3<double, size_t>(&bus, 3);

  sys->addComponent(comp1);
  sys->addComponent(comp3);
  sys->addComponent(bus1);
  sys->addNode(&bus);

  sys->allocate();

  std::vector<double> y  = {1, 2, 3, 6, 7, 8, 4, 5};
  std::vector<double> yp = {1, 2, 3, 6, 7, 8, 4, 5};

  auto sys_y  = sys->y().getData();
  auto sys_yp = sys->yp().getData();

  for (size_t i = 0; i < sys->size(); i++)
  {
    sys_y[i]  = y[i];
    sys_yp[i] = yp[i];
  }

  sys->evaluateResidual();

  GridKit::HiresComponent3<double, size_t>*       comp3copy    = new GridKit::HiresComponent3<double, size_t>(*comp3);
  GridKit::BusPartitionInterface<double, size_t>* busInterface = new GridKit::BusPartitionInterface<double, size_t>(&bus, comp3copy, 4);

  partition1->addComponent(comp1);
  partition1->addComponent(bus1);
  partition1->addInterface(busInterface);
  partition1->addNode(&bus);

  partition2->addComponent(comp3);

  auto subsystems = {partition1, partition2};

  for (auto* partition : subsystems)
  {
    partition->allocate();
    partition->initialize();
    partition->updateTime(2, 5);

    for (size_t i = 0; i < partition->getExternSize(); i++)
    {
      partition->getExternalDataY()[i]  = y[partition->getExternalDataIndices()[i]];
      partition->getExternalDataYP()[i] = yp[partition->getExternalDataIndices()[i]];
    }

    auto* partition_y  = partition->y().getData();
    auto* partition_yp = partition->yp().getData();

    for (size_t i = 0; i < partition->getInternalSize(); i++)
    {
      partition_y[i]  = y[partition->getNodeConnection(i)];
      partition_yp[i] = yp[partition->getNodeConnection(i)];
    }

    partition->y().setDataUpdated();
    partition->yp().setDataUpdated();

    partition->evaluateResidual();
  }

  printPartitionResiduals(partition1, "Partition 1");
  printPartitionResiduals(partition2, "Partition 2");

  printHeader("Reference Output");
  auto* residual = sys->getResidual().getData();

  for (size_t i = 0; i < sys->getInternalSize(); ++i)
  {
    printRow(residual[i], i);
  }

  delete sys;
  delete partition1;
  delete partition2;

  return 0;
}

void printRow(double residual, size_t index)
{
  std::cout << std::left << std::setw(10)
            << std::defaultfloat << std::setprecision(5)
            << residual
            << "  ----------  "
            << std::right << std::setw(10)
            << index
            << '\n';
}

void printHeader(const std::string& title)
{
  std::cout << "\n------------- " << title << " -----------\n";

  std::cout << std::left << std::setw(10) << "Res Values"
            << "  ----------  "
            << std::right << std::setw(10) << "Comp. Index"
            << '\n';
}

void printPartitionResiduals(auto* partition, const std::string& title)
{
  printHeader(title);

  auto* residual = partition->getResidual().getData();

  for (size_t i = 0; i < partition->getInternalSize(); ++i)
  {
    const auto comp_index = partition->getNodeConnection(i);
    printRow(residual[i], comp_index);
  }
}