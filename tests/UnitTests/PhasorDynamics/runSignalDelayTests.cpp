#include "SignalDelayTests.hpp"

int main()
{
  GridKit::Testing::TestingResults   results;
  GridKit::Testing::SignalDelayTests tests;
  results += tests.prehistoryDerivative();
  results += tests.derivativeCorner();
  results += tests.interpolation();
  results += tests.historyLifecycle();
  results += tests.compression();
  results += tests.continuityAndHiddenFeature();
  results += tests.component();
  results += tests.declarations();
  return results.summary();
}
