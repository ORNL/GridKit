#include "ConverterRegcaTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::ConverterRegcaTests<double, size_t> test;

  result += test.constructionAndValidation();
  result += test.parameterValidation();
  result += test.initializesFromPowerFlowAndPublishesSignals();
  result += test.unconnectedCommandsRemainConstant();
  result += test.initializesAboveHighVoltageLimit();
  result += test.rejectsInitializationBelowLvacmBreakpoint();
  result += test.rejectsInitializationWithActiveLvacm();
  result += test.initializesAtLvacmUpperBreakpoint();
  result += test.rejectsZeroTerminalVoltage();
  result += test.signalVerification();
  result += test.nullBusVerification();
  result += test.busInjectionUsesSystemBase();
  result += test.residualEquations();
  result += test.highVoltageReactiveCurrentRoot();
  result += test.positiveInitialReactivePowerSelectsUpperIqRateLimit();
  result += test.disabledLvplRemovesIlDependence();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif

  return result.summary();
}
