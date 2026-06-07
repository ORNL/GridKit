#include "DelayTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::DelayTests<double, size_t> test;

  result += test.historyInterpolates();
  result += test.historyPrunesToWindow();
  result += test.delayShiftsInput();
  result += test.delayPrehistoryBeforeStart();
  result += test.delayConstantPrehistory();
  result += test.delayMaxStepSize();
  result += test.delayRejectsNonpositiveDelay();

  return result.summary();
}
