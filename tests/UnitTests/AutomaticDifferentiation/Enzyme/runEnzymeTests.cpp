/**
 * @file runEnzymeTests.cpp
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 *
 */

#include "EnzymeTests.hpp"

int main()
{
  using namespace GridKit;
  using namespace GridKit::Testing;

  GridKit::Testing::TestingResults              result;
  GridKit::Testing::EnzymeTests<double, size_t> test;

  result += test.scalar_square();

  return result.summary();
}
