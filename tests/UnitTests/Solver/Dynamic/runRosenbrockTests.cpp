#include "RosenbrockTests.hpp"

int main()
{
  using namespace GridKit;
  using namespace GridKit::Testing;

  GridKit::Testing::TestingResults               result;
  GridKit::Testing::RosenbrockTests<double, int> test;

  result += test.test();

  return result.summary();
}
