#include <cstddef>

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>
#include <iostream>

#include <GridKit/Model/PowerElectronics/PartitionInterface/BusPartitionInterface.hpp>
#include <GridKit/Model/PowerElectronics/PartitionInterface/ComponentPartitionInterface.hpp>
#include <GridKit/Model/PowerElectronics/SystemModelPowerElectronics.hpp>

#include "HiresBus.hpp"
#include "HiresComponent1.hpp"
#include "HiresComponent3.hpp"

std::vector<double> hires(const std::vector<double>& y, const std::vector<double>& yp);

int main(int /* argc */, char const** /* argv */)
{

  double abs_tol         = 1.0e-5;
  double rel_tol         = 1.0e-5;
  size_t max_step_number = 3000;
  bool   use_jac         = true;

  auto* partition1 = new GridKit::PowerElectronicsModel<double, size_t>(rel_tol, abs_tol, use_jac, max_step_number);
  auto* partition2 = new GridKit::PowerElectronicsModel<double, size_t>(rel_tol, abs_tol, use_jac, max_step_number);

  GridKit::HiresComponent1<double, size_t>* comp1 = new GridKit::HiresComponent1<double, size_t>(1);
  GridKit::HiresBus<double, size_t>*        bus1  = new GridKit::HiresBus<double, size_t>(2);
  GridKit::HiresComponent3<double, size_t>* comp3 = new GridKit::HiresComponent3<double, size_t>(3);

  comp1->setExternalConnectionNodes(0, 0);
  comp1->setExternalConnectionNodes(1, 1);
  comp1->setExternalConnectionNodes(2, 2);
  comp1->setExternalConnectionNodes(3, 3);
  comp1->setExternalConnectionNodes(4, 4);

  bus1->setExternalConnectionNodes(0, 3);
  bus1->setExternalConnectionNodes(1, 4);

  comp3->setExternalConnectionNodes(0, 5);
  comp3->setExternalConnectionNodes(1, 6);
  comp3->setExternalConnectionNodes(2, 7);
  comp3->setExternalConnectionNodes(3, 3);
  comp3->setExternalConnectionNodes(4, 4);

  GridKit::HiresComponent3<double, size_t>              comp2copy(*comp3);
  GridKit::BusPartitionInterface<double, size_t>*       busInterface  = new GridKit::BusPartitionInterface<double, size_t>(comp2copy, 3, 4, 4);
  GridKit::ComponentPartitionInterface<double, size_t>* compInterface = new GridKit::ComponentPartitionInterface<double, size_t>(comp3, 3, 4, 5);

  partition1->addComponent(comp1);
  partition1->addComponent(bus1);
  partition1->addComponent(busInterface);

  partition2->addComponent(compInterface);

  size_t sys_size = 8;

  partition1->allocate(sys_size);
  partition2->allocate(sys_size);

  busInterface->setExternalDataY({6, 7, 8});
  busInterface->setExternalDataYP({6, 7, 8});

  compInterface->setExternalDataY({4, 5});
  compInterface->setExternalDataYP({4, 5});

  for (size_t i = 0; i < sys_size; i++)
  {

    partition1->y()[i]  = static_cast<double>(i) + 1;
    partition1->yp()[i] = static_cast<double>(i) + 1;

    partition2->y()[i]  = static_cast<double>(i) + 1;
    partition2->yp()[i] = static_cast<double>(i) + 1;
  }

  partition1->evaluateResidual();
  partition2->evaluateResidual();

  auto ref = hires({1, 2, 3, 4, 5, 6, 7, 8}, {1, 2, 3, 4, 5, 6, 7, 8});

  std::cout << "------------- Partition 1 -----------" << std::endl;
  for (size_t i = 0; i < sys_size; i++)
  {
    printf("%-10.5g  ----------  %10.5g\n", partition1->getResidual()[i], ref[i]);
  }

  std::cout << "\n------------- Partition 2 -----------" << std::endl;
  for (size_t i = 0; i < sys_size; i++)
  {
    printf("%-10.5g  ----------  %10.5g\n", partition2->getResidual()[i], ref[i]);
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
