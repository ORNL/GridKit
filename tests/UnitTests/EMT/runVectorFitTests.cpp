#include "VectorFitTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::VectorFitTests<double, size_t> test;

  result += test.validation();
  result += test.residual();
  result += test.steadyState();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif

  return result.summary();
}
