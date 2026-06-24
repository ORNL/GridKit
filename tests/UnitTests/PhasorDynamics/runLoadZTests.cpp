#include "LoadZTests.hpp"

int main()
{
  using namespace GridKit;
  using namespace GridKit::Testing;

  GridKit::Testing::TestingResults             result;
  GridKit::Testing::LoadZTests<double, size_t> test;

  result += test.constructor();
  result += test.residual();
  result += test.jacobian();
  result += test.monitor();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.enzymeJacobian();
#endif

  return result.summary();
}
