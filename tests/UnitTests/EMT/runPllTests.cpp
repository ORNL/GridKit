#include "PllTests.hpp"

int main()
{
  GridKit::Testing::TestingResults           result;
  GridKit::Testing::PllTests<double, size_t> test;
  result += test.wiring();
  result += test.initialState();
  result += test.residual();
  result += test.lockIn();
  result += test.phaseJump();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif
  return result.summary();
}
