#include "LoadZIPTests.hpp"

int main()
{
  using namespace GridKit;
  using namespace GridKit::Testing;

  GridKit::Testing::TestingResults               result;
  GridKit::Testing::LoadZIPTests<double, size_t> test;

  result += test.constructor();
  result += test.residual();

  return result.summary();
}
