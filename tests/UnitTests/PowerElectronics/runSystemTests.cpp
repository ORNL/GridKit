#include "SystemTests.hpp"

int main()
{
  using namespace GridKit;
  using namespace GridKit::Testing;

  GridKit::Testing::TestingResults              result;
  GridKit::Testing::SystemTests<double, size_t> test;

  result += test.jacobianSparsity();

  return result.summary();
}
