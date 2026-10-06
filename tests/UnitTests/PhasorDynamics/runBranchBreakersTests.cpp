#include "BranchBreakersTests.hpp"

int main()
{
  GridKit::Testing::TestingResults                      result;
  GridKit::Testing::BranchBreakersTests<double, size_t> test;

  result += test.residual();
  result += test.latch();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
  result += test.infiniteBusJacobian();
#endif

  return result.summary();
}
