#include "SystemTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::SystemTests<double, size_t> test;

  result += test.parse();
  result += test.system();
  result += test.boundaryAssembly();
#ifdef GRIDKIT_ENABLE_SUNDIALS
  result += test.steadyState();
  result += test.threeBusSteadyState();
  result += test.recursiveSteadyState();
  result += test.switchEnergization();
  result += test.machineFlatStart();
  result += test.machineGovernorFlatStart();
  result += test.twinCircuit();
  result += test.twinLine();
  result += test.twinSource();
#endif

  return result.summary();
}
