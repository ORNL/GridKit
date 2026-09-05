#include "Ieeet1Tests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;
  GridKit::Testing::Ieeet1Tests    test;
  result += test.initialization();
  result += test.saturation();
  result += test.validation();
  result += test.residualAndLimits();
  result += test.jacobian();
  result += test.repeatedJacobian();
  result += test.dependencyTracking();
  return result.summary();
}
