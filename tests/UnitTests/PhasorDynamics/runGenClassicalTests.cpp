#include "GenClassicalTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::GenClassicalTests<double, size_t> test;

  result += test.constructor();
  result += test.residual();
  result += test.initial();
  result += test.zeroInitialResidual();

  return result.summary();
}
