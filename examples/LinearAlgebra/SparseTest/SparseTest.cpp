

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <tuple>
#include <vector>

#include <GridKit/LinearAlgebra/SparseMatrix/COO_Matrix.hpp>

int main()
{
  std::vector<double> val{0.1, 0.2, 0.3, 0.4};
  std::vector<size_t> x{2, 1, 3, 1};
  std::vector<size_t> y{1, 3, 2, 2};
  size_t              n = 4;
  size_t              m = 4;

  GridKit::LinearAlgebra::COO_Matrix<double, size_t> A = GridKit::LinearAlgebra::COO_Matrix<double, size_t>(x, y, val, m, n);

  std::vector<double> valn(4);
  std::vector<size_t> xn(4);
  std::vector<size_t> yn(4);

  std::tie(xn, yn, valn) = A.getEntries();

  for (size_t i = 0; i < valn.size(); i++)
  {
    std::cout << valn[i] << "\n";
  }

  std::cout << "A:\n";
  A.printMatrix();

  std::vector<double>                                val2{0.5, 0.6, 0.7, 0.8, 1.0};
  std::vector<size_t>                                x2{0, 2, 0, 2, 1};
  std::vector<size_t>                                y2{3, 3, 2, 2, 3};
  GridKit::LinearAlgebra::COO_Matrix<double, size_t> B = GridKit::LinearAlgebra::COO_Matrix<double, size_t>(x2, y2, val2, m, n);

  std::cout << "B:\n";
  B.printMatrix();

  A.axpy(2.0, B);

  std::cout << "A + 2B:\n";
  A.printMatrix();

  std::vector<size_t> r;
  std::vector<size_t> c;
  std::vector<double> v;
  std::tie(r, c, v) = A.setDataToCSR();

  for (size_t i = 0; i < r.size() - 1; i++)
  {
    std::cout << r[i] << std::endl;
    size_t rdiff = r[i + 1] - r[i];
    for (size_t j = 0; j < rdiff; j++)
    {
      std::cout << c[j + r[i]] << ", " << v[j + r[i]] << std::endl;
    }
  }
  std::cout << r[r.size() - 1] << std::endl;

  // Basic Verification test
  std::vector<size_t> rtest   = {0, 2, 4, 7, 8};
  std::vector<size_t> ctest   = {2, 3, 2, 3, 1, 2, 3, 2};
  std::vector<double> valtest = {1.4, 1.0, 0.4, 2.2, 0.1, 1.6, 1.2, 0.3};

  assert(rtest.size() == r.size());
  assert(ctest.size() == c.size());
  assert(valtest.size() == v.size());

  int failval = 0;
  for (size_t i = 0; i < rtest.size(); i++)
  {
    if (r[i] != rtest[i])
    {
      failval--;
    }
  }
  for (size_t i = 0; i < ctest.size(); i++)
  {
    double vdiff = v[i] - valtest[i];
    if (c[i] != ctest[i] || -1e-14 > vdiff || vdiff > 1e-14)
    {
      failval--;
    }
  }

  if (failval == 0)
  {
    std::cout << "Success!" << std::endl;
  }
  else
  {
    std::cout << "Failed!" << std::endl;
  }

  std::vector<double> v3 = {1.2, 1.8, 1.0, 1.6, 1.4, 1.7, 1.1, 1.5, 1.3};
  std::vector<size_t> x3 = {1, 2, 0, 0, 2, 1, 0, 2, 1};
  std::vector<size_t> y3 = {1, 2, 0, 0, 0, 1, 1, 2, 2};
  size_t              m3 = 3;
  size_t              n3 = 3;

  GridKit::LinearAlgebra::COO_Matrix<double, size_t> C;
  C.resetEntries(x3, y3, v3, m3, n3);

  std::vector<size_t> csr_r;
  std::vector<size_t> csr_c;
  std::vector<double> csr_v;
  std::tie(csr_r, csr_c, csr_v) = C.getCsrData();

  std::cout << "C after getCsrData:\n";
  C.printMatrix();

  for (size_t i = 0; i < csr_r.size() - 1; i++)
  {
    std::cout << csr_r[i] << std::endl;
    size_t rdiff = csr_r[i + 1] - csr_r[i];
    for (size_t j = 0; j < rdiff; j++)
    {
      std::cout << csr_c[j + csr_r[i]] << ", " << csr_v[j + csr_r[i]] << std::endl;
    }
  }
  std::cout << csr_r[csr_r.size() - 1] << std::endl;

  // Basic Verification test
  std::vector<size_t> csr_rtest   = {0, 2, 4, 6};
  std::vector<size_t> csr_ctest   = {0, 1, 1, 2, 0, 2};
  std::vector<double> csr_valtest = {2.6, 1.1, 2.9, 1.3, 1.4, 3.3};
  std::vector<size_t> maptest     = {2, 5, 0, 0, 4, 2, 1, 5, 3};

  std::vector<size_t>& map_to_csr = C.getMapToCsr();

  assert(csr_rtest.size() == csr_r.size());
  assert(csr_ctest.size() == csr_c.size());
  assert(csr_valtest.size() == csr_v.size());
  assert(maptest.size() == map_to_csr.size());

  for (size_t i = 0; i < csr_rtest.size(); i++)
  {
    if (csr_r[i] != csr_rtest[i])
    {
      failval--;
    }
  }
  for (size_t i = 0; i < csr_ctest.size(); i++)
  {
    double vdiff = csr_v[i] - csr_valtest[i];
    if (csr_c[i] != csr_ctest[i] || -1e-14 > vdiff || vdiff > 1e-14)
    {
      failval--;
    }
  }

  for (size_t i = 0; i < maptest.size(); i++)
  {
    if (map_to_csr[i] != maptest[i])
    {
      failval--;
    }
  }

  if (failval == 0)
  {
    std::cout << "Success!" << std::endl;
  }
  else
  {
    std::cout << "Failed!" << std::endl;
  }

  return failval;
}
