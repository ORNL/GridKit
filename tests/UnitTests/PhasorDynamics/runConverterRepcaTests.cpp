#include "ConverterRepcaTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::ConverterRepcaTests<double, size_t> test;

  result += test.validation();
  result += test.initializationAndSignals();
  result += test.initializationDomain();
  result += test.residualEquations();
  result += test.reactiveControl();
  result += test.activePowerControl();
  result += test.derivatives();
  result += test.jacobian();

  return result.summary();
}
