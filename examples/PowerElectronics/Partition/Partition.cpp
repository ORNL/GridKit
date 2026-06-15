#include <cstddef>

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>

#include <GridKit/Model/PowerElectronics/PartitionInterface/BusPartitionInterface.hpp>
#include <GridKit/Model/PowerElectronics/SubsystemModel.hpp>

#include "HiresBus.hpp"
#include "HiresComponent1.hpp"
#include "HiresComponent3.hpp"

std::vector<double> hires(const std::vector<double>& y, const std::vector<double>& yp);

int main(int /* argc */, char const** /* argv */)
{

  GridKit::SubsystemModel<double, size_t>* partition1 = new GridKit::SubsystemModel<double, size_t>();
  GridKit::SubsystemModel<double, size_t>* partition2 = new GridKit::SubsystemModel<double, size_t>();

  GridKit::HiresComponent1<double, size_t>* comp1 = new GridKit::HiresComponent1<double, size_t>(1);
  GridKit::HiresBus<double, size_t>*        bus1  = new GridKit::HiresBus<double, size_t>(2);
  GridKit::HiresComponent3<double, size_t>* comp3 = new GridKit::HiresComponent3<double, size_t>(3);

  comp1->setExternalConnectionNodes(0, 3);
  comp1->setExternalConnectionNodes(1, 4);
  comp1->setExternalConnectionNodes(2, 0);
  comp1->setExternalConnectionNodes(3, 1);
  comp1->setExternalConnectionNodes(4, 2);

  bus1->setExternalConnectionNodes(0, 3);
  bus1->setExternalConnectionNodes(1, 4);

  comp3->setExternalConnectionNodes(0, 3);
  comp3->setExternalConnectionNodes(1, 4);
  comp3->setExternalConnectionNodes(2, 5);
  comp3->setExternalConnectionNodes(3, 6);
  comp3->setExternalConnectionNodes(4, 7);

  GridKit::HiresComponent3<double, size_t>        comp3copy(*comp3);
  GridKit::BusPartitionInterface<double, size_t>* busInterface = new GridKit::BusPartitionInterface<double, size_t>(comp3copy, 3, 4, 4);

  partition1->addComponent(comp1);
  partition1->addComponent(bus1);
  partition1->addComponent(busInterface);

  partition2->addComponent(comp3);

  partition1->allocate();
  partition2->allocate();

  std::vector<double> y  = {1, 2, 3, 4, 5, 6, 7, 8};
  std::vector<double> yp = {1, 2, 3, 4, 5, 6, 7, 8};

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
    partition1->y()[i]  = static_cast<double>(1 + i);
    partition1->yp()[i] = static_cast<double>(1 + i);
  }

  // Distribute internals to partition 2
  for (size_t i = 0; i < partition2->getInternalSize(); i++)
  {
    partition2->y()[i]  = static_cast<double>(6 + i);
    partition2->yp()[i] = static_cast<double>(6 + i);
  }

  // Evaluate Residuals for each partition
  partition1->evaluateResidual();
  partition2->evaluateResidual();

  auto printTitle = [](std::string msg) -> void
  {
    std::cout << "\n------------- " << msg << " -----------" << std::endl;
    printf("%-10s  ----------  %10s\n", "Res Values", "Comp. Index");
  };

  // Print Residuals from partition 1
  printTitle("Partition 1");
  for (size_t i = 0; i < partition1->getInternalSize(); i++)
  {
    auto com_index = partition1->getNodeConnection(static_cast<size_t>(i));
    printf("%-10.5g  ----------  %10zu\n", partition1->getResidual()[i], com_index);
  }

  // Print Residuals from partition 2
  printTitle("Partition 2");
  for (size_t i = 0; i < partition2->getInternalSize(); i++)
  {
    auto com_index = partition2->getNodeConnection(static_cast<size_t>(i));
    printf("%-10.5g  ----------  %10zu\n", partition2->getResidual()[i], com_index);
  }

  auto ref = hires(y, yp);

  printTitle("Reference Solution");
  for (size_t i = 0; i < 8; i++)
  {
    printf("%-10.5g  ----------  %10zu\n", ref[i], i);
  }

  delete partition1;
  delete partition2;
}

std::vector<double> hires(const std::vector<double>& y, const std::vector<double>& yp)
{
  std::vector<double> f(8);

  f[0] = yp[0] + 1.71 * y[0] - 0.43 * y[1] - 8.32 * y[2] - 0.0007;
  f[1] = yp[1] - 1.71 * y[0] + 8.75 * y[1];
  f[2] = yp[2] + 10.03 * y[2] - 0.43 * y[3] - 0.035 * y[4];
  f[3] = yp[3] - 8.32 * y[1] - 1.71 * y[2] + 1.12 * y[3];
  f[4] = yp[4] + 1.745 * y[4] - 0.43 * y[5] - 0.43 * y[6];
  f[5] = yp[5] + 280 * y[5] * y[7] - 0.69 * y[3] - 1.71 * y[4] + 0.43 * y[5] - 0.69 * y[6];
  f[6] = yp[6] - 280 * y[5] * y[7] + 1.81 * y[6];
  f[7] = yp[7] + 280 * y[5] * y[7] - 1.81 * y[6];

  return f;
}
