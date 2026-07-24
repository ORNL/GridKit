// #include <cstddef>

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <string>

#include <GridKit/Model/PowerElectronics/Bus/MicrogridBus.hpp>
#include <GridKit/Model/PowerElectronics/PartitionInterface/BusPartitionInterface.hpp>
#include <GridKit/Model/PowerElectronics/SubsystemModel.hpp>

#include "HiresBus.hpp"
#include "HiresComponent1.hpp"
#include "HiresComponent3.hpp"

std::vector<double> hires(const std::vector<double>& y, const std::vector<double>& yp);

int main(int /* argc */, char const** /* argv */)
{

  using Bus = GridKit::PowerElectronics::MicrogridBus<double, size_t>;
  Bus bus;

  GridKit::SubsystemModel<double, size_t>* partition1 = new GridKit::SubsystemModel<double, size_t>();
  GridKit::SubsystemModel<double, size_t>* partition2 = new GridKit::SubsystemModel<double, size_t>();

  GridKit::HiresComponent1<double, size_t>* comp1 = new GridKit::HiresComponent1<double, size_t>(&bus, 1);
  GridKit::HiresBus<double, size_t>*        bus1  = new GridKit::HiresBus<double, size_t>(&bus, 2);
  GridKit::HiresComponent3<double, size_t>* comp3 = new GridKit::HiresComponent3<double, size_t>(&bus, 3);

  bus.allocate();

  bus.setExternalConnectionNodes(0, 3);
  bus.setExternalConnectionNodes(1, 4);

  comp1->allocate();
  comp3->allocate();
  bus1->allocate();

  comp1->setExternalConnectionNodes(2, 0);
  comp1->setExternalConnectionNodes(3, 1);
  comp1->setExternalConnectionNodes(4, 2);

  bus1->setExternalConnectionNodes(0, 3);
  bus1->setExternalConnectionNodes(1, 4);

  comp3->setExternalConnectionNodes(2, 5);
  comp3->setExternalConnectionNodes(3, 6);
  comp3->setExternalConnectionNodes(4, 7);

  GridKit::HiresComponent3<double, size_t>        comp3copy(*comp3);
  GridKit::BusPartitionInterface<double, size_t>* busInterface = new GridKit::BusPartitionInterface<double, size_t>(bus, comp3copy, 4);

  busInterface->allocate();
  partition1->addComponent(comp1);
  partition1->addComponent(bus1);
  partition1->addComponent(busInterface);
  partition1->addNode(&bus);

  partition2->addComponent(comp3);

  auto subsystems = {partition1, partition2};

  std::vector<double> y  = {1, 2, 3, 4, 5, 6, 7, 8};
  std::vector<double> yp = {1, 2, 3, 4, 5, 6, 7, 8};

  // Distribute externals to partition 1

  for (auto* partition : subsystems)
  {
    partition->allocate();
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
    partition->evaluateJacobian();
  }

  auto printHeader = [](const std::string& title)
  {
    std::cout << "\n------------- " << title << " -----------\n";

    std::cout << std::left << std::setw(10) << "Res Values"
              << "  ----------  "
              << std::right << std::setw(10) << "Comp. Index"
              << '\n';
  };

  auto printRow = [](double residual, size_t index)
  {
    std::cout << std::left << std::setw(10)
              << std::defaultfloat << std::setprecision(5)
              << residual
              << "  ----------  "
              << std::right << std::setw(10)
              << index
              << '\n';
  };

  auto printPartitionResiduals = [&](auto* partition, const std::string& title)
  {
    printHeader(title);

    auto* residual = partition->getResidual().getData();

    for (size_t i = 0; i < partition->getInternalSize(); ++i)
    {
      const auto comp_index = partition->getNodeConnection(i);
      printRow(residual[i], comp_index);
    }
  };

  printPartitionResiduals(partition1, "Partition 1");
  printPartitionResiduals(partition2, "Partition 2");

  auto ref = hires(y, yp);

  printHeader("Reference Output");
  for (size_t i = 0; i < ref.size(); ++i)
  {
    printRow(ref[i], i);
  }

  delete comp1;
  delete bus1;
  delete busInterface;
  delete comp3;

  delete partition1;
  delete partition2;

  return 0;
}

std::vector<double> hires(const std::vector<double>& y, const std::vector<double>& yp)
{
  std::vector<double> f(8);

  f[0] = -yp[0] - 1.71 * y[0] + 0.43 * y[1] + 8.32 * y[2] + 0.0007;
  f[1] = -yp[1] + 1.71 * y[0] - 8.75 * y[1];
  f[2] = -yp[2] - 10.03 * y[2] + 0.43 * y[3] + 0.035 * y[4];
  f[3] = -yp[3] + 8.32 * y[1] + 1.71 * y[2] - 1.12 * y[3];
  f[4] = -yp[4] - 1.745 * y[4] + 0.43 * y[5] + 0.43 * y[6];
  f[5] = -yp[5] - 280 * y[5] * y[7] + 0.69 * y[3] + 1.71 * y[4] - 0.43 * y[5] + 0.69 * y[6];
  f[6] = -yp[6] + 280 * y[5] * y[7] - 1.81 * y[6];
  f[7] = -yp[7] - 280 * y[5] * y[7] + 1.81 * y[6];

  return f;
}

// template <typename RealT, typename IdxT>
// int jac_hires(const std::vector<double>& y, [[maybe_unused]] const std::vector<double>& yp, GridKit::LinearAlgebra::DenseMatrix<RealT, IdxT>& jac, RealT alpha)
// {

//   jac.setValue(0, 0, -alpha - 1.71);
//   jac.setValue(0, 1, 0.43);
//   jac.setValue(0, 2, 8.32);

//   jac.setValue(1, 0, 1.71);
//   jac.setValue(1, 1, -alpha - 8.75);

//   jac.setValue(2, 2, -alpha - 10.03);
//   jac.setValue(2, 3, 0.43);
//   jac.setValue(2, 4, 0.035);

//   jac.setValue(3, 1, 8.32);
//   jac.setValue(3, 2, 1.71);
//   jac.setValue(3, 3, -alpha - 1.12);

//   jac.setValue(4, 4, -alpha - 1.745);
//   jac.setValue(4, 5, 0.43);
//   jac.setValue(4, 6, 0.43);

//   jac.setValue(5, 3, 0.69);
//   jac.setValue(5, 4, 1.71);
//   jac.setValue(5, 5, -alpha - 280.0 * y[7] - 0.43);
//   jac.setValue(5, 6, 0.69);
//   jac.setValue(5, 7, -280.0 * y[5]);

//   jac.setValue(6, 5, 280.0 * y[7]);
//   jac.setValue(6, 6, -alpha - 1.81);
//   jac.setValue(6, 7, 280.0 * y[5]);

//   jac.setValue(7, 5, -280.0 * y[7]);
//   jac.setValue(7, 6, 1.81);
//   jac.setValue(7, 7, -alpha - 280.0 * y[5]);

//   return 0;
// }
