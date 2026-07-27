#include "GovernorHygovTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::GovernorHygovTests<double, size_t> test;

  result += test.validation();
  result += test.initializationAndSignals();
  result += test.initializationDomain();
  result += test.residualEquations();
  result += test.governorControl();
  result += test.turbineDynamics();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif

  return result.summary();
}
