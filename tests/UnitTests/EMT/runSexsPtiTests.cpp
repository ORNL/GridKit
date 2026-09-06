#include "SexsPtiTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;
  GridKit::Testing::SexsPtiTests   test;
  result += test.initializationAndMeasurement();
  result += test.phasorParityAndLimits();
  result += test.jacobianAndDependencies();
  result += test.computedSignalsAndPattern();
  result += test.validation();
  return result.summary();
}
