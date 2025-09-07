#include "CsrTests.hpp"

int main()
{
  using namespace GridKit;
  using namespace GridKit::Testing;
  using namespace GridKit::LinearAlgebra;

  TestingResults             result;
  CsrTests<double, unsigned> test;

  result += test.cooToCsrTest();
  result += test.testCooMove();
  result += test.testCsrBuilderTemplate();
  result += test.testCsrBuilderComplete();
  result += test.testUnsortedMatrix();
  result += test.testUnsortedToSorted();

  return result.summary();
}
