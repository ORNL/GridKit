#include "SystemTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;
  GridKit::Testing::SystemTests    test;

  result += test.allocation();
  result += test.reference();
  result += test.jacobian();
  result += test.gradient();
  result += test.hessian();
  result += test.solve();
  result += test.solutionState();
  result += test.parseMatpowerData();
  result += test.applyMatpowerData();

  return result.summary();
}
