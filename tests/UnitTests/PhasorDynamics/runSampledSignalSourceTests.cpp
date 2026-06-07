#include "SampledSignalSourceTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::SampledSignalSourceTests<double, size_t> test;

  result += test.signalHistory();
  result += test.senderPrehistory();
  result += test.inlineSource();
  result += test.csvSource();

  return result.summary();
}
