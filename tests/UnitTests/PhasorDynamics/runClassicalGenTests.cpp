#include "ClassicalGenTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::ClassicalGenTests<double, size_t> test;

  result += test.constructor();
  result += test.residual();

  return result.summary();
}