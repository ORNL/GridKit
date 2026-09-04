#include "MachineTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::MachineTests<double, size_t> test;

  result += test.initialization();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif

  return result.summary();
}
