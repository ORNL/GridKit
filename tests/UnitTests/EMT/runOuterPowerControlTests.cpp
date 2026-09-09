#include "OuterPowerControlTests.hpp"

int main()
{
  GridKit::Testing::TestingResults                         result;
  GridKit::Testing::OuterPowerControlTests<double, size_t> test;
  result += test.wiring();
  result += test.initialState();
  result += test.residual();
  result += test.powerMeasurements();
  result += test.steadyState();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif
  return result.summary();
}
