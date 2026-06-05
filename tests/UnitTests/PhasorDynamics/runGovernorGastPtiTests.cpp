#include "GovernorGastPtiTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::GovernorGastPtiTests<double, size_t> test;

  result += test.constructor();
  result += test.zeroInitialResidual();
  result += test.prefSignal();
  result += test.residual();
  result += test.antiWindupLimiter();
  result += test.timeConstantTags();
  result += test.parameterValidation();
  result += test.signalValidation();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif

  return result.summary();
}
