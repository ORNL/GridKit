#include "SignalDelayTests.hpp"

int main()
{
  GridKit::Testing::TestingResults   results;
  GridKit::Testing::SignalDelayTests tests;
  results += tests.interpolation();
  results += tests.historyLifecycle();
  results += tests.compression();
  results += tests.component();
  results += tests.declarations();
  return results.summary();
}
