#include "SystemTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::SystemTests<double, size_t> test;

  result += test.parse();
  result += test.system();
#ifdef GRIDKIT_ENABLE_SUNDIALS
  result += test.steadyState();
  result += test.threeBusSteadyState();
  result += test.switchEnergization();
  result += test.machineFlatStart();
  result += test.machineGovernorFlatStart();
#endif

  return result.summary();
}
