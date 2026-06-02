/**
 * @file runForcedOscillationTests.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Test runner for the forced oscillation signal operator.
 */

#include "ForcedOscillationTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::ForcedOscillationTests<double, size_t> test;

  result += test.constructor();
  result += test.source_mode();
  result += test.additive_mode();
  result += test.delayed_additive_identity();
  result += test.consistent_initial_condition();
  result += test.sine_chirp_decay();
  result += test.raised_cosine_window();
  result += test.clamp_defaults_and_limits();
  result += test.monitor_values();
  result += test.validation();
  result += test.json_and_system_wiring();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif

  return result.summary();
}
