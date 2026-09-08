#include "TransformerTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::TransformerTests<double, size_t> test;

  result += test.wiring();
  result += test.residual();
  result += test.steadyState();
  result += test.connection();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif

  return result.summary();
}
