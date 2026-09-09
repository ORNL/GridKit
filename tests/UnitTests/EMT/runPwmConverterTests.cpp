#include "PwmConverterTests.hpp"

int main()
{
  GridKit::Testing::TestingResults    result;
  GridKit::Testing::PwmConverterTests test;
  result += test.waveform();
  result += test.validation();
  result += test.bridgeVoltages();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.signalGradients();
#endif
  result += test.dependencyTracking();
  result += test.parseAndAssemble();
  result += test.constantSignals();
  result += test.runtimeSmoothing();
  result += test.continuousMean();
  result += test.continuousResolution();
  result += test.continuousInput();
  result += test.voltageCommand();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.commandGradients();
#endif
  result += test.signalReadScope();
  result += test.monitors();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
  result += test.computedPwmJacobian();
#endif
#ifdef GRIDKIT_ENABLE_SUNDIALS
  result += test.integrate();
#endif
  return result.summary();
}
