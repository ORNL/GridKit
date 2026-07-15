#include "GovernorGastPtiTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::GovernorGastPtiTests<double, size_t> test;

  result += test.constructor();
  result += test.zeroInitialResidual();
  result += test.baseConversion();
  result += test.absoluteTolerance();
  result += test.prefSignal();
  result += test.prefSignalBaseConversion();
  result += test.residual();
  result += test.antiWindupLimiter();
  result += test.responseModes();
  result += test.narrowLimitWarning();
  result += test.initializationValidation();
  result += test.smoothMinimumInitialization();
  result += test.smoothMinimumEqualityRejected();
  result += test.timeConstantMinimum();
  result += test.parameterValidation();
  result += test.signalValidation();
  result += test.jsonParseAndSystemAssembly();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif

  return result.summary();
}
