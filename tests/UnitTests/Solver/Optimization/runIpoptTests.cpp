#include "IpoptTests.hpp"

int main()
{
  using namespace GridKit::Testing;
  using namespace GridKit::Testing::Optimization;

  TestingResults                  result;
  IpoptTests<double, std::size_t> test;

  result += test.exactIpoptCallbacks();
  result += test.emptyExactHessian();
  result += test.secondOrderDerivativeTest();

  return result.summary();
}
