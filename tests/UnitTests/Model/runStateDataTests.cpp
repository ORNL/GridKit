#include "StateDataTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;
  GridKit::Testing::StateDataTests test;

  result += test.keys();
  result += test.roundTrip();
  result += test.rejectsText();
  result += test.terminalConversion();

  return result.summary();
}
