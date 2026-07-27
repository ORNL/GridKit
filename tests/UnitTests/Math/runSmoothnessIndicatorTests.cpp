#include "SmoothnessIndicatorTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::SmoothnessIndicatorTests<double> test;

  result += test.clamp();
  result += test.deadband1();
  result += test.deadband2();
  result += test.limitIndicators();
  result += test.slew();
  result += test.linseg();
  result += test.ramp();
  result += test.minMax();
  result += test.antiWindupIndicator();
  result += test.antiWindup();
  result += test.dynamicAntiWindupBounds();

  return result.summary();
}
