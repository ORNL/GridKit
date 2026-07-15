#include "ConverterReecbTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::ConverterReecbTests<double, size_t> test;

  result += test.validation();
  result += test.signals();
  result += test.publishRefs();
  result += test.baseSignals();
  result += test.feedbackBase();
  result += test.zeroTime();
  result += test.qPriority();
  result += test.pPriority();
  result += test.voltageBand();
  result += test.piSaturation();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif
  result += test.json();

  return result.summary();
}
