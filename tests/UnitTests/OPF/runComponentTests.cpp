#include "BranchComponentTests.hpp"
#include "ComponentTests.hpp"

int main()
{
  GridKit::Testing::TestingResults          result;
  GridKit::Testing::OPFComponentTests       component_tests;
  GridKit::Testing::OPFBranchComponentTests branch_tests;

  result += component_tests.busAndLoadExactEmptyStructures();
  result += component_tests.generatorExactDerivatives();
  result += component_tests.shuntExactDerivatives();
  result += branch_tests.terminalPowerAndDerivativeStructures();
  result += branch_tests.exactJacobianAndHessian();

  return result.summary();
}
