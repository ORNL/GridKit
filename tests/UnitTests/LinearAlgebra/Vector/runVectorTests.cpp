#include <fstream>
#include <iostream>
#include <string>

#include "VectorTests.hpp"

int main(int, char**)
{
  GridKit::Testing::TestingResults result;

  {
    GridKit::Testing::VectorTests<double, size_t> test;

    std::cout << "Running vector tests on CPU:\n";
    result += test.vectorConstructor(50, 5);
    result += test.vectorConstructor(50);

    result += test.setData(50);

    result += test.copyFromExternal(50);
    // result += test.copyToExternal(50);

    result += test.resize(100, 50);
    result += test.resizeEmpty(5);

    result += test.setToConst(50);
    result += test.setToConst(50);
  }

#ifdef GRIDKIT_ENABLE_GPU
  {
    GridKit::Testing::VectorTests test(GridKit::memory::DEVICE);

    std::cout << "Running Testing on GPU:\n";
    result += test.vectorConstructor(50, 5);
    result += test.vectorConstructor(50);

    result += test.setData(50);

    result += test.copyFromExternal(50);
    // result += test.copyToExternal(50);

    result += test.resize(100, 50);

    result += test.setToConst(50);
    result += test.setToConst(50);
    result += test.syncData(50);
    result += test.syncData(50);
  }
#endif

  return result.summary();
}
