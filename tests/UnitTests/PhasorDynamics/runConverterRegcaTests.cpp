#include "ConverterRegcaTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::ConverterRegcaTests<double, size_t> test;

  result += test.constructionAndValidation();
  result += test.parameterValidation();
  result += test.initializesFromPowerFlowAndPublishesSignals();
  result += test.baseSignals();
  result += test.unconnectedCommandsRemainConstant();
  result += test.externalCommandsDriveRuntimeResidual();
  result += test.initializesAboveHighVoltageLimit();
  result += test.initializesBelowLvacmBreakpoint();
  result += test.initializesOnLvacmRampWithActivePower();
  result += test.rejectsZeroTerminalVoltage();
  result += test.signalVerification();
  result += test.nullBusVerification();
  result += test.busInjectionUsesSystemBase();
  result += test.residualEquations();
  result += test.highVoltageReactiveCurrentRoot();
  result += test.limiterBranchCoverage();
  result += test.jsonParseAndSystemAssembly();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif

  return result.summary();
}
