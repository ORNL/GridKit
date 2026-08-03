#include "StabilizerIeeestTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::StabilizerIeeestTests<double, size_t> test;

  result += test.constructor();
  result += test.zeroInitialResidual();
  result += test.residual();
  result += test.verifyNegativeParameters();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif

  return result.summary();
}
