#include "ModulationTests.hpp"

int main()
{
  GridKit::Testing::TestingResults  result;
  GridKit::Testing::ModulationTests test;
  result += test.limits();
  result += test.gradients();
  return result.summary();
}
