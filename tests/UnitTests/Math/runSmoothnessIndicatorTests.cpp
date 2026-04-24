#include "SmoothnessIndicatorTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::SmoothnessIndicatorTests<double> test;

  result += test.antiWindupIndicator();

  return result.summary();
}
