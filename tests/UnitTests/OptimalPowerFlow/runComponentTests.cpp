#include "ComponentTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;
  GridKit::Testing::ComponentTests test;

  result += test.branchPower();
  result += test.branchJacobian();
  result += test.branchHessian();
  result += test.branchPattern();
  result += test.generator();
  result += test.offlineGenerator();
  result += test.shunt();
  result += test.bus();
  result += test.load();

  return result.summary();
}
