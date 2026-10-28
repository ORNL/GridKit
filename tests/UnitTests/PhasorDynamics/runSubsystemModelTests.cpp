#include "SubsystemModelTests.hpp"

int main()
{
  GridKit::Testing::TestingResults                      result;
  GridKit::Testing::SubsystemModelTests<double, size_t> test;

  result += test.flatVsHierarchical();
  result += test.hierarchicalInput();
  result += test.invalidHierarchy();
  result += test.systemBases();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif

  return result.summary();
}
