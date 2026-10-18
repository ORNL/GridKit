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
  result += test.outOfServiceGenerator();
  result += test.bus();
  result += test.busShunt();
  result += test.loadConstantPower();
  result += test.loadEqualities();
  result += test.loadJacobian();
  result += test.loadHessian();
  result += test.loadZeroVoltage();
  result += test.loadConfiguration();

  return result.summary();
}
