#include "SystemComputedOutputTests.hpp"

int main()
{
  using namespace GridKit;
  using namespace GridKit::Testing;

  GridKit::Testing::TestingResults                            result;
  GridKit::Testing::SystemComputedOutputTests<double, size_t> test;

  result += test.outputsPrecedeResiduals();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.chainRuleJacobian();
#endif

  return result.summary();
}
