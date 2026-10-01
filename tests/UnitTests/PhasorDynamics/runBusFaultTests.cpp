#include "BusFaultTests.hpp"

int main()
{
  using namespace GridKit;
  using namespace GridKit::Testing;

  GridKit::Testing::TestingResults                result;
  GridKit::Testing::BusFaultTests<double, size_t> test;

  result += test.constructor();
  result += test.zeroInitialResidual(true);
  result += test.zeroInitialResidual(false);
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian(true);
  result += test.jacobian(false);
#endif

  return result.summary();
}
