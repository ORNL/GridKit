#include "GasCompressorTests.hpp"

int main()
{
  GridKit::Testing::TestingResults                     result;
  GridKit::Testing::GasCompressorTests<double, size_t> test;

  result += test.constructor();
  result += test.parameterValidation();
  result += test.initialization();
  result += test.residual();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
  result += test.fixedPressureBus();
#endif

  return result.summary();
}
