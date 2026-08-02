#include "GovernorGastPtiTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::GovernorGastPtiTests<double, size_t> test;

  result += test.validation();
  result += test.initializationAndSignals();
  result += test.initializationDomain();
  result += test.initializationExactness();
  result += test.residualEquations();
  result += test.governorControl();
  result += test.temperatureLimiting();
  result += test.jacobian();

  return result.summary();
}
