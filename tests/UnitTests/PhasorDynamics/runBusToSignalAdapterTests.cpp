#include "BusToSignalAdapterTests.hpp"

int main()
{
  using namespace GridKit;
  using namespace GridKit::Testing;

  TestingResults                          result;
  BusToSignalAdapterTests<double, size_t> test;

  result += test.constructor();
  result += test.residual();

  return result.summary();
}
