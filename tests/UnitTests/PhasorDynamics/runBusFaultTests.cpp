#include "BusFaultTests.hpp"

int main()
{
  using namespace GridKit;
  using namespace GridKit::Testing;

  GridKit::Testing::TestingResults                result;
  GridKit::Testing::BusFaultTests<double, size_t> test;

  result += test.constructor();
  result += test.zeroInitialResidual(true);
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian(true);
  result += test.jacobianAfterClearing();
#endif

  return result.summary();
}
