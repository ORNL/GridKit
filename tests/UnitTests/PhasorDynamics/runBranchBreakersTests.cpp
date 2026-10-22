#include "BranchBreakersTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::BranchBreakersTests<double, size_t> test;

  result += test.validation();
  result += test.closedMatchesBranch();
  result += test.openEnds();
  result += test.latchResidual();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif

  return result.summary();
}
