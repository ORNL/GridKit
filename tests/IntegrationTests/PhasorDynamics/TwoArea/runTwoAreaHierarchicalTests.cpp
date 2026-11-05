/**
 * @file runTwoAreaHierarchicalTests.cpp
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 * @brief Driver for the TwoArea hierarchical system test.
 *
 */

#include "TwoAreaHierarchicalTests.hpp"

int main()
{
  GridKit::Testing::TestingResults                           result;
  GridKit::Testing::TwoAreaHierarchicalTests<double, size_t> test;

  result += test.flatVsHierarchical();

  return result.summary();
}
