#include "CSRTests.hpp"

int main()
{
  using namespace GridKit;
  using namespace GridKit::Testing;
  using namespace GridKit::LinearAlgebra;

  TestingResults             result;
  CSRTests<double, unsigned> test;

  result += test.cooToCsrTest();

  return result.summary();
}