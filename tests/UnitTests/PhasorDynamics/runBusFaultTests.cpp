#include "BusFaultTests.hpp"

int main()
{
  using namespace GridKit;
  using namespace GridKit::Testing;

  GridKit::Testing::TestingResults                result;
  GridKit::Testing::BusFaultTests<double, size_t> test;

  result += test.constructor();
  result += test.residual(true);
  result += test.residual(false);
  result += test.jacobian(true);
  result += test.jacobian(false);
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.enzymeJacobian(true);
  result += test.enzymeJacobian(false);
  result += test.jacobianAcrossStatusChanges();
#endif

  return result.summary();
}
