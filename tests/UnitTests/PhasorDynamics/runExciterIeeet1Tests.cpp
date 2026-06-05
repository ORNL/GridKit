#include "ExciterIeeet1Tests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::ExciterIeeet1Tests<double, size_t> test;

  result += test.constructor();
  result += test.zeroInitialResidual();
  result += test.zeroTrIsAlgebraic();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif

  return result.summary();
}
