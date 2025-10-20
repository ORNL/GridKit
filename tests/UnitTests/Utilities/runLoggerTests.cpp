/**
 * @file runLoggerTests.cpp
 * @brief Driver for Logger class tests.
 * @author Slaven Peles <peless@ornl.org>
 */

#include <iostream>
#include <ostream>

#include "LoggerTests.hpp"

int main()
{
  using namespace GridKit;

  // Create LoggerTests object
  GridKit::Testing::LoggerTests test;

  // Create test results accounting object
  GridKit::Testing::TestingResults result;

  // Run tests
  result += test.errorOutput();
  result += test.warningOutput();
  result += test.summaryOutput();
  result += test.miscOutput();

  // Return tests summary
  return result.summary();
}
