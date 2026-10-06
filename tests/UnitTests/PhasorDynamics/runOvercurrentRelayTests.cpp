#include "OvercurrentRelayTests.hpp"

int main()
{
  GridKit::Testing::TestingResults                        result;
  GridKit::Testing::OvercurrentRelayTests<double, size_t> test;

  result += test.residual();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif

  return result.summary();
}
