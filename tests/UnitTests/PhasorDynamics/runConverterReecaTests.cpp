#include "ConverterReecaTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::ConverterReecaTests<double, size_t> test;

  result += test.constructor();
  result += test.parameterValidation();
  result += test.initialization();
  result += test.residual();
  result += test.zeroTimeConstantTags();
  result += test.outputAvailability();
  result += test.signalVerification();
  result += test.deviceInitFallbacks();
  result += test.priorityInitialization();
  result += test.unsupportedInitializationRejects();
  result += test.jsonParseAndSystemAssembly();
  result += test.initFallbackSystemAssembly();
  result += test.systemRejectsUnlinkedRequiredSignals();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif

  return result.summary();
}
