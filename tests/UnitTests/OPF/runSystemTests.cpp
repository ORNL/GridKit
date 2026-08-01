#include "SystemTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;
  GridKit::Testing::OPFSystemTests test;

  result += test.rejectInvalidSystems();
  result += test.allocateInitializeAndBounds();
  result += test.rejectNonfiniteEvaluations();
  result += test.evaluateAndDifferentiate();
  result += test.writeSolutionState();

  return result.summary();
}
