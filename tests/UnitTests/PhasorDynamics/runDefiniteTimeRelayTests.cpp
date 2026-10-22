#include "DefiniteTimeRelayTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::DefiniteTimeRelayTests<double, size_t> test;

  result += test.validation();
  result += test.residualEquations();
  result += test.initialization();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
  result += test.breakerCurrentJacobian();
#endif

  return result.summary();
}
