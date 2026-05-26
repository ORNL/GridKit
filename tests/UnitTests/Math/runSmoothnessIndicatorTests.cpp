#include "SmoothnessIndicatorTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::SmoothnessIndicatorTests<double> test;

  result += test.clamp();
  result += test.deadband();
  result += test.limitIndicators();
  result += test.slew();
  result += test.linseg();
  result += test.ramp();
  result += test.minMax();
  result += test.antiWindupIndicator();
  result += test.antiWindup();

  return result.summary();
}
