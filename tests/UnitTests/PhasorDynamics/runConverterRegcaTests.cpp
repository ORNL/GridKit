#include "ConverterRegcaTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::ConverterRegcaTests<double, size_t> test;

  result += test.constructor();
  result += test.parameterValidation();
  result += test.lifecycle();
  result += test.steadyStateInitializationGolden();
  result += test.attachedSignalInitialization();
  result += test.invalidInitialization();
  result += test.signalVerification();
  result += test.nullBusVerification();
  result += test.residualGoldenVectors();
  result += test.busInjection();
  result += test.jsonParseAndSystemAssembly();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif

  return result.summary();
}
