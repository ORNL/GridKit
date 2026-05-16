#include "GenrouTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::GenrouTests<double, size_t> test;

  result += test.constructor();
  result += test.accessors();
  result += test.hard_coded_residual();
  result += test.residual();
  result += test.monitor_system_base();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif

  return result.summary();
}
