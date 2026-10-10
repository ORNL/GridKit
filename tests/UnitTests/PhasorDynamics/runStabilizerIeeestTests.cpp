#include "StabilizerIeeestTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::StabilizerIeeestTests<double, size_t> test;

  result += test.validation();
  result += test.initializationAndSignals();
  result += test.residualEquations();
  result += test.monitor();

#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif

  return result.summary();
}
