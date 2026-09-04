#include "DependentVoltageSourceTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::DependentVoltageSourceTests<double, size_t> test;

  result += test.seriesResidual();
  result += test.rationalResidual();
  result += test.validation();
  result += test.parseAndBuild();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif

  return result.summary();
}
