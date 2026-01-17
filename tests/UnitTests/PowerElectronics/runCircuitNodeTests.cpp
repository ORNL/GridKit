#include "CircuitNodeTests.hpp"

int main()
{
  GridKit::Testing::TestingResults            result;
  GridKit::Testing::NodeTests<double, size_t> test;

  result += test.constructor();
  result += test.residual();

  return result.summary();
}