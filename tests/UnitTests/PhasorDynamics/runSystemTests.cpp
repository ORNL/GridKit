#include "SystemTests.hpp"

int main()
{
  using namespace GridKit;
  using namespace GridKit::Testing;

  GridKit::Testing::TestingResults              result;
  GridKit::Testing::SystemTests<double, size_t> test;

  result += test.constructor();
  result += test.composer();
  result += test.reallocateAfterTopologyChange();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif

  result += test.allocationError();
  result += test.signalError();

  return result.summary();
}
