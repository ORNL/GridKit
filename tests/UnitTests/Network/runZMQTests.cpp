#include "ZMQTests.hpp"

int main()
{
  using namespace GridKit::Testing;

  TestingResults           result;
  ZMQTests<double, size_t> test;

  result += test.basic();

  return result.summary();
}
