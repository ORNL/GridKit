#include "SystemTests.hpp"

int main()
{
  using namespace GridKit;
  using namespace GridKit::Testing;

  GridKit::Testing::TestingResults              result;
  GridKit::Testing::SystemTests<double, size_t> test;

  result += test.jacobianSparsity(2);
  result += test.jacobianSparsity(128);

  return result.summary();
}
