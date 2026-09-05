#include "PwmConverterTests.hpp"

int main()
{
  GridKit::Testing::TestingResults    result;
  GridKit::Testing::PwmConverterTests test;
  result += test.waveform();
  result += test.validation();
  result += test.bridgeVoltages();
  result += test.signalGradients();
  result += test.dependencyTracking();
  result += test.parseAndAssemble();
  result += test.monitors();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif
#ifdef GRIDKIT_ENABLE_SUNDIALS
  result += test.integrate();
#endif
  return result.summary();
}
