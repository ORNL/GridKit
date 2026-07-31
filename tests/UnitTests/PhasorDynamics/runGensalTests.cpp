#include "GensalTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::GensalTests<double, size_t> test;

  result += test.constructor();
  result += test.hard_coded_residual();
  result += test.residual();
  result += test.residual_nonzero_ra();
  result += test.frequency_base();
  result += test.monitor_system_base();
  result += test.operating_state_signals();

#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif
  return result.summary();
}
