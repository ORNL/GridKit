#include "GastPtiTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;
  GridKit::Testing::GastPtiTests   test;
  result += test.initializationAndBases();
  result += test.phasorParityAndTemperatureLimit();
  result += test.jacobianAndComputedSignals();
  result += test.validation();
  return result.summary();
}
