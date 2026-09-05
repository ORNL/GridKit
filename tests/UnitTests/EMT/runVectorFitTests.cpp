#include "VectorFitTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::VectorFitTests<double, size_t> test;

  result += test.validation();
  result += test.residual();
  result += test.steadyState();
  result += test.jacobian();

  return result.summary();
}
