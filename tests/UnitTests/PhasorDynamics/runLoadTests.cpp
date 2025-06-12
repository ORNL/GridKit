#include "LoadTests.hpp"

int main()
{
  using namespace GridKit;
  using namespace GridKit::Testing;

  GridKit::Testing::TestingResults            result;
  GridKit::Testing::LoadTests<double, size_t> test;

  result += test.constructor();
  result += test.residual();
  result += test.jacobian();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.enzyme_jacobian();
#endif

  return result.summary();
}
