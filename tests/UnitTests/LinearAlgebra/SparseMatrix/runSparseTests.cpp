#include "SparseTests.hpp"
#include <LinearAlgebra/LinAlgWorkspaceCpu.hpp>

using namespace GridKit;
using namespace LinearAlgebra;
using namespace Testing;

/**
 * @brief Run sparse matrix tests with a given backend
 *
 * @tparam WorkspaceType workspace type LinAlgWorkspace{Cpu, CUDA, HIP} supported
 * @param[in] backend - name of the hardware backend
 * @param[out] result - test results
 */
template <typename WorkspaceType>
void runTests(const std::string& backend, memory::MemorySpace memspace, TestingResults& result)
{
  std::cout << "Running tests on " << backend << ":\n";

  WorkspaceType workspace;
  workspace.initializeHandles();

  SparseTests test(memspace);

  //   ReSolve::io::Logger::setVerbosity(ReSolve::io::Logger::NONE);

  result += test.constructor(50, 50, 2);
  result += test.constructor(50, 100, 2);
  result += test.constructor(50, 50, 2, true);
  result += test.constructor(50, 50, 2, false, false);

  result += test.setDataPointers(50, 50, 2);
  result += test.setDataPointers(50, 100, 2);
  result += test.setValuesPointer(50, 50, 2);
  result += test.setValuesPointer(50, 100, 2);

  result += test.copyValues(50, 50, 2);
  result += test.copyValues(50, 100, 2);
  result += test.copyValuesAndSetValues(50, 50, 2);
  result += test.copyValuesAndSetValues(50, 100, 2);
  result += test.copyValuesAndSetDataPointers(50, 50, 2);
  result += test.copyValuesAndSetDataPointers(50, 100, 2);

  result += test.allocateAndDestroyData(50, 50, 2);
  result += test.allocateAndDestroyData(50, 100, 2);
}

int main(int, char**)
{
  TestingResults result;
  runTests<LinAlgWorkspaceCpu>("CPU", memory::HOST, result);

  return result.summary();
}
