#include "CSRTests.hpp"

int main()
{
  using namespace GridKit;
  using namespace GridKit::Testing;
  using namespace GridKit::LinearAlgebra;

  TestingResults             result;
  CSRTests<double, unsigned> test;

  result += test.cooToCsrTest();
  result += test.testCooMove();
  result += test.testCsrBuilderTemplate();
  result += test.testCsrBuilderComplete();

  return result.summary();
}