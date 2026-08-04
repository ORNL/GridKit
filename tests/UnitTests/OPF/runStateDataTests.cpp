#include "StateDataTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;
  GridKit::Testing::StateDataTests test;

  result += test.parseCompleteState();
  result += test.parseOmittedAndNullValues();
  result += test.parseLegacyGoldenState();
  result += test.rejectInvalidInput();
  result += test.writeDeterministicState();
  result += test.writeCompleteState();
  result += test.fileAndOutputFailures();
  result += test.rejectNonfiniteOutput();

  return result.summary();
}
