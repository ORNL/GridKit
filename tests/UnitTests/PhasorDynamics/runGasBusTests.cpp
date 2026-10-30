#include "GasBusTests.hpp"

int main()
{
  GridKit::Testing::TestingResults              result;
  GridKit::Testing::GasBusTests<double, size_t> test;

  result += test.constructor();
  result += test.initialization();
  result += test.residual();
  result += test.withdrawal();
  result += test.parameterValidation();
  result += test.fixedPressure();

  return result.summary();
}
