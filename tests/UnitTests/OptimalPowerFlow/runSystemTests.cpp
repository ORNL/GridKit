#include "SystemTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;
  GridKit::Testing::SystemTests    test;

  result += test.allocation();
  result += test.jacobian();
  result += test.gradient();
  result += test.hessianLinearity();
  result += test.unboundedInequalities();
  result += test.configuration();
  result += test.parseCase();
  result += test.connections();
  result += test.objectiveAssembly();
  result += test.equalityAssembly();
  result += test.inequalityAssembly();
  result += test.hessianAssembly();
  result += test.repeatedVariables();
  result += test.connectionOwnership();

  return result.summary();
}
