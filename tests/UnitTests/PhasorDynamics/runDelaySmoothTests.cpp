#include "DelaySmoothTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::DelaySmoothTests<double, size_t> test;

  result += test.derivedBlockCount();
  result += test.verifyRejectsBadParams();
  result += test.initializeFlattensChain();
  result += test.residualZeroAtSteadyState();
  result += test.residualRecursion();
  result += test.jacobianBanded();
  result += test.noMaxStepCap();

  return result.summary();
}
