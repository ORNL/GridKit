#include "SparseCooTests.hpp"

using namespace GridKit;
using namespace LinearAlgebra;
using namespace Testing;

/**
 * @brief Run COO sparse matrix tests with a given backend
 *
 * @param[in] backend  - name of the hardware backend
 * @param[in] memspace - memory space for the tests
 * @param[out] result  - test results
 */
template <class ScalarT, typename IdxT>
void runTests(const std::string& backend, memory::MemorySpace memspace, TestingResults& result)
{
  std::cout << "Running tests on " << backend << ":\n";

  SparseCooTests<ScalarT, IdxT> test(memspace);

  result += test.constructor(50, 50, 2);
  result += test.constructor(50, 100, 2);

  result += test.setDataPointers(50, 50, 2);
  result += test.setDataPointers(50, 100, 2);

  result += test.getCsrRowData();
}

int main(int, char**)
{
  TestingResults result;
  runTests<double, size_t>("CPU", memory::HOST, result);
  runTests<double, long int>("CPU", memory::HOST, result);

  return result.summary();
}
