#include "LoadShuntBranchTests.hpp"

int main()
{
  GridKit::Testing::TestingResults       result;
  GridKit::Testing::LoadShuntBranchTests test;

  result += test.loadBehavior();
  result += test.shuntBehavior();
  result += test.shuntFiniteDifferenceJacobian();
  result += test.branchBehavior();
  result += test.branchFiniteDifferenceJacobian();

  return result.summary();
}
