#include "LoadZIPTests.hpp"

int main()
{
  using namespace GridKit;
  using namespace GridKit::Testing;

  GridKit::Testing::TestingResults               result;
  GridKit::Testing::LoadZIPTests<double, size_t> test;

  result += test.constructor();
  result += test.initialization();
  result += test.signalInputs();
  result += test.reconnectableJacobianStructure();
  result += test.residual();
  result += test.monitor();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif

  return result.summary();
}
