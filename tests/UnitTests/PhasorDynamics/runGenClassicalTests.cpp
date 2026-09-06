#include "GenClassicalTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::GenClassicalTests<double, size_t> test;

  result += test.constructor();
  result += test.initial();
  result += test.residual();
  result += test.residual_nonzero_ra();
  result += test.frequency_base();
  result += test.monitor_system_base();
  result += test.signals();
  result += test.hard_coded_residual();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif

  return result.summary();
}
