#include "LineLumpedTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::LineLumpedTests<double, size_t> test;

  result += test.wiring();
  result += test.residual();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif

  return result.summary();
}
