#include "InverseTimeRelayTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::InverseTimeRelayTests<double, size_t> test;

  result += test.validation();
  result += test.residualEquations();
  result += test.initialization();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif

  return result.summary();
}
