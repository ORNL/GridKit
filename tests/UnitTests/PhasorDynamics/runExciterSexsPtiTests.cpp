#include "ExciterSexsPtiTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::ExciterSexsPtiTests<double, size_t> test;

  result += test.constructor();
  result += test.zeroInitialResidual();
  result += test.vsSignal();
  result += test.antiWindupLimiter();
  result += test.parameterValidation();
  result += test.systemAssembly();
  result += test.systemSignalFreshness();

  return result.summary();
}
