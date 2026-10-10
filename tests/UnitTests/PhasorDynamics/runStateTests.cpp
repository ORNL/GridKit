#include "StateTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;
  GridKit::Testing::StateTests     test;

  result += test.applyState();
  result += test.removals();
  result += test.offlineMachine();

  return result.summary();
}
