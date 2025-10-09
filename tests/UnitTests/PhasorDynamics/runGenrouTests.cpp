#include "GenrouTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::GenrouTests<double, size_t> test;

  result += test.constructor();
  result += test.accessors();
  result += test.hard_coded_residual();
  result += test.residual();

#if 0 // Disabled GenrouEnzyme because of an intermittent build issue
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif
#endif
  return result.summary();
}
