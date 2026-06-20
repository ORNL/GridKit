#include "ArkodeTests.hpp"

int main()
{
  GridKit::Testing::TestingResults              result;
  GridKit::Testing::ArkodeTests<double, size_t> test;

  result += test.arkStepIdentityMass();
  result += test.arkStepNonIdentityMass();
  result += test.erkStepForwardEuler();
  result += test.erkStepRejectsNonIdentityMass();
  result += test.callback();

  return result.summary();
}
