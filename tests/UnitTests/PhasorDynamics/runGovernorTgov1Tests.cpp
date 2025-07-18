#include "GovernorTgov1Tests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::GovernorTgov1Tests<double, size_t> test;

  result += test.constructor();
  result += test.accessors();
  result += test.residual();

  return result.summary();
}
