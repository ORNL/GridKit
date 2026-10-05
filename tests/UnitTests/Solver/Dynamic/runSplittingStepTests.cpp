#include "SplittingStepTests.hpp"

int main()
{
  GridKit::Testing::TestingResults     result;
  GridKit::Testing::SplittingStepTests test;
  result += test.sequentialRamps();
  result += test.event();
  result += test.rejectionAndThreads();
  result += test.regionalFailure();
  result += test.configuration();
  result += test.diagnostics();
  return result.summary();
}
