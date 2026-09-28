#include "SourceDependentNortonTests.hpp"

int main()
{
  using namespace GridKit;
  using namespace GridKit::Testing;

  GridKit::Testing::TestingResults                             result;
  GridKit::Testing::SourceDependentNortonTests<double, size_t> test;

  result += test.validation();
  result += test.initialization();
  result += test.residual();
  result += test.monitor();
  result += test.jacobian();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.enzymeJacobian();
#endif

  return result.summary();
}
