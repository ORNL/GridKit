#include <fstream>
#include <iostream>
#include <string>

#include "VectorHandlerTests.hpp"

int main(int, char**)
{
  GridKit::Testing::TestingResults result;

  {
    std::cout << "Running vector handler tests on CPU:\n";

    GridKit::LinearAlgebra::VectorHandler<double, size_t> handler;

    GridKit::Testing::VectorHandlerTests<double, size_t> test(handler);
    result += test.vectorHandlerConstructor();
    result += test.dot(50);
    result += test.axpy(50);
    result += test.scal(50);
    result += test.amax(50);
    result += test.gemv(500, 10);
    result += test.axpyMulti(100, 10);
    result += test.massDot(100, 10);
    result += test.scale(100);
    result += test.diagSolve(100);
    result += test.max(100);
    result += test.abs(100);

    std::cout << "\n";
  }

#ifdef GRIDKIT_ENABLE_CUDA
  // TODO: Once a CUDA VectorHandler implementation is available, run the
  // same tests here with a device-backed VectorHandler and memory::DEVICE.
#endif

#ifdef GRIDKIT_ENABLE_HIP
  // TODO: Once a HIP VectorHandler implementation is available, run the
  // same tests here with a device-backed VectorHandler and memory::DEVICE.
#endif

  return result.summary();
}
