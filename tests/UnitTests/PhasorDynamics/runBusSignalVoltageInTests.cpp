#include "BusSignalVoltageInTests.hpp"

int main()
{
  using namespace GridKit;
  using namespace GridKit::Testing;

  TestingResults                          result;
  BusSignalVoltageInTests<double, size_t> test;

  result += test.constructor();
  result += test.voltageInputs();
  result += test.currentOutputs();
  result += test.verifyUnlinked();
  result += test.dependencyTracking();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif

  return result.summary();
}
