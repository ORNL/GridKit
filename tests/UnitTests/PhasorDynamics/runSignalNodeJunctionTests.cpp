/**
 * @file runSignalNodeJunctionTests.cpp
 * @brief Driver for signal-node junction unit tests.
 */
#include <cstddef>

#include "SignalNodeJunctionTests.hpp"

int main()
{
  GridKit::Testing::TestingResults                               result;
  GridKit::Testing::SignalNodeJunctionTests<double, std::size_t> test;

  result += test.algebraicEquation();
  result += test.designatedInputInitialization();
  result += test.invalidConfiguration();
  result += test.initializationCycle();

  return result.summary();
}
