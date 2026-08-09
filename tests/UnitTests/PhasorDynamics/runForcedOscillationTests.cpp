/**
 * @file runForcedOscillationTests.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Test runner for the forced-oscillation signal source.
 */

#include <GridKit/Utilities/Logger/Logger.hpp>

#include "ForcedOscillationTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  using Log                     = GridKit::Utilities::Logger;
  const auto previous_verbosity = Log::verbosity();
  Log::setVerbosity(Log::Verbosity::NONE);

  GridKit::Testing::ForcedOscillationTests<double, size_t> test;
  result += test.validation();
  result += test.initialization();
  result += test.waveform();
  result += test.carrierWaveforms();
  result += test.activationWindowAndMonitors();
  result += test.dependencyTracking();

  Log::setVerbosity(previous_verbosity);
  return result.summary();
}
