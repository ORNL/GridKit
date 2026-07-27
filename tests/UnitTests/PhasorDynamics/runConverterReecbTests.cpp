#include "ConverterReecbTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::ConverterReecbTests<double, size_t> test;

  result += test.validation();
  result += test.initializationAndSignals();
  result += test.initializationDomain();
  result += test.residualEquations();
  result += test.reactiveControl();
  result += test.activePowerControl();
  result += test.currentPriority();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif

  return result.summary();
}
