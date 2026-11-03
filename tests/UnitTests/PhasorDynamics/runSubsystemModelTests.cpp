/**
 * @file runSubsystemModelTests.cpp
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 * @brief Driver for hierarchical system model tests.
 *
 */

#include "SubsystemModelTests.hpp"

int main()
{
  GridKit::Testing::TestingResults                      result;
  GridKit::Testing::SubsystemModelTests<double, size_t> test;

  result += test.flatVsHierarchical();
  result += test.hierarchicalInput();
  result += test.hierarchicalLookup();
  result += test.emptyIdLookup();
  result += test.duplicateIds();
  result += test.invalidHierarchy();
  result += test.systemBases();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif

  return result.summary();
}
