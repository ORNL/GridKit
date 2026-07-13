#include <cstddef>

#include "GridKit/LinearAlgebra/DenseMatrix/DenseMatrix.hpp"

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

template <typename RealT, typename IdxT>
int jac_hires(const std::vector<double>& y, [[maybe_unused]] const std::vector<double>& yp, GridKit::LinearAlgebra::DenseMatrix<RealT, IdxT>& jac, RealT alpha);

template <typename RealT, typename IdxT>
void addCsrToDense(
    GridKit::LinearAlgebra::CsrMatrix<RealT, IdxT>*   csr,
    GridKit::LinearAlgebra::DenseMatrix<RealT, IdxT>& dense,
    GridKit::SubsystemModel<double, size_t>*          partition);

int main(int /* argc */, char const** /* argv */)
{

  using Bus = GridKit::PowerElectronics::MicrogridBus<double, size_t>;
  Bus bus;

  GridKit::SubsystemModel<double, size_t>* partition1 = new GridKit::SubsystemModel<double, size_t>();
  GridKit::SubsystemModel<double, size_t>* partition2 = new GridKit::SubsystemModel<double, size_t>();

  GridKit::HiresComponent1<double, size_t>* comp1 = new GridKit::HiresComponent1<double, size_t>(1);
  GridKit::HiresBus<double, size_t>*        bus1  = new GridKit::HiresBus<double, size_t>(2);
  GridKit::HiresComponent3<double, size_t>* comp3 = new GridKit::HiresComponent3<double, size_t>(3);

  bus.allocate();

  bus.setExternalConnectionNodes(0, 3);
  bus.setExternalConnectionNodes(1, 4);

  comp1->allocate();
  comp3->allocate();
  bus1->allocate();

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

  GridKit::HiresComponent3<double, size_t>
                                                  comp3copy(*comp3);
  GridKit::BusPartitionInterface<double, size_t>* busInterface = new GridKit::BusPartitionInterface<double, size_t>(bus, comp3copy, 4);

  busInterface->allocate();
  partition1->addComponent(comp1);
  partition1->addComponent(bus1);
  partition1->addComponent(busInterface);
  partition1->addNode(&bus);

  partition2->addComponent(comp3);

  std::cout << "Msg" << std::endl;
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

  partition1->updateTime(0, 1.0);
  partition2->updateTime(0, 1.0);

  partition1->evaluateJacobian();
  partition2->evaluateJacobian();

  GridKit::LinearAlgebra::DenseMatrix<double, size_t> jac1(8, 8);
  GridKit::LinearAlgebra::DenseMatrix<double, size_t> jac2(8, 8);

  GridKit::LinearAlgebra::DenseMatrix<double, size_t> jac_ref(8, 8);

  addCsrToDense(partition1->getCsrJacobian(), jac1, partition1);
  addCsrToDense(partition2->getCsrJacobian(), jac2, partition2);
  jac_hires(y, yp, jac_ref, 1.0);

  auto print_dense = [&](GridKit::LinearAlgebra::DenseMatrix<double, size_t>& jc) -> void
  {
    for (size_t i = 0; i < 8; ++i)
    {
      for (size_t j = 0; j < 8; ++j)
      {
        std::cout << std::setw(10) << jc.getValue(i, j) << ' ';
      }
      std::cout << '\n';
    }
    std::cout << '\n';
  };

  print_dense(jac1);
  print_dense(jac2);
  print_dense(jac_ref);

  delete comp1;
  delete bus1;
  delete busInterface;
  delete comp3;
  delete partition1;
  delete partition2;
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

template <typename RealT, typename IdxT>
int jac_hires(const std::vector<double>& y, [[maybe_unused]] const std::vector<double>& yp, GridKit::LinearAlgebra::DenseMatrix<RealT, IdxT>& jac, RealT alpha)
{

  jac.setValue(0, 0, -alpha - 1.71);
  jac.setValue(0, 1, 0.43);
  jac.setValue(0, 2, 8.32);

  jac.setValue(1, 0, 1.71);
  jac.setValue(1, 1, -alpha - 8.75);

  jac.setValue(2, 2, -alpha - 10.03);
  jac.setValue(2, 3, 0.43);
  jac.setValue(2, 4, 0.035);

  jac.setValue(3, 1, 8.32);
  jac.setValue(3, 2, 1.71);
  jac.setValue(3, 3, -alpha - 1.12);

  jac.setValue(4, 4, -alpha - 1.745);
  jac.setValue(4, 5, 0.43);
  jac.setValue(4, 6, 0.43);

  jac.setValue(5, 3, 0.69);
  jac.setValue(5, 4, 1.71);
  jac.setValue(5, 5, -alpha - 280.0 * y[7] - 0.43);
  jac.setValue(5, 6, 0.69);
  jac.setValue(5, 7, -280.0 * y[5]);

  jac.setValue(6, 5, 280.0 * y[7]);
  jac.setValue(6, 6, -alpha - 1.81);
  jac.setValue(6, 7, 280.0 * y[5]);

  jac.setValue(7, 5, -280.0 * y[7]);
  jac.setValue(7, 6, 1.81);
  jac.setValue(7, 7, -alpha - 280.0 * y[5]);

  return 0;
}

template <typename RealT, typename IdxT>
void addCsrToDense(
    GridKit::LinearAlgebra::CsrMatrix<RealT, IdxT>*   csr,
    GridKit::LinearAlgebra::DenseMatrix<RealT, IdxT>& dense,
    GridKit::SubsystemModel<double, size_t>*          partition)
{
  auto* row_ptr = csr->getRowData(GridKit::LinearAlgebra::memory::HOST);
  auto* col_ind = csr->getColData(GridKit::LinearAlgebra::memory::HOST);
  auto* values  = csr->getValues(GridKit::LinearAlgebra::memory::HOST);

  assert(row_ptr != nullptr);
  assert(col_ind != nullptr);
  assert(values != nullptr);

  for (IdxT row = 0; row < csr->getNumRows(); ++row)
  {
    for (IdxT k = row_ptr[row]; k < row_ptr[row + 1]; ++k)
    {
      const IdxT col = col_ind[k];

      auto r = partition->getNodeConnection(row);
      auto c = partition->getNodeConnection(col);

      dense.setValue(r, c, dense.getValue(r, c) + values[k]);
    }
  }
}
