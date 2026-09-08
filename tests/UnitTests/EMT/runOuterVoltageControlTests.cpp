#include "OuterVoltageControlTests.hpp"

int main()
{
  GridKit::Testing::TestingResults                           result;
  GridKit::Testing::OuterVoltageControlTests<double, size_t> test;
  result += test.outerControl();
  return result.summary();
}
