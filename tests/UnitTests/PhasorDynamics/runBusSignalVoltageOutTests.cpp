#include "BusSignalVoltageOutTests.hpp"

int main()
{
  using namespace GridKit;
  using namespace GridKit::Testing;

  TestingResults                           result;
  BusSignalVoltageOutTests<double, size_t> test;

  result += test.constructor();
  result += test.voltageOutputs();
  result += test.storageBinding();
  result += test.residual();
  result += test.verifyUnlinked();
  result += test.dependencyTracking();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif

  return result.summary();
}
