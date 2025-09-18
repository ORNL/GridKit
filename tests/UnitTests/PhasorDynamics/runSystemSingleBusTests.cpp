#include "SystemSingleBusTests.hpp"

int main()
{
  using namespace GridKit;
  using namespace GridKit::Testing;

  GridKit::Testing::TestingResults                 result;
  GridKit::Testing::SystemSingleBusTests<double, size_t> test;

  result += test.constructor();

  return result.summary();
}
