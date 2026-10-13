#include "SignalDelayIntegrationTests.hpp"

int main()
{
  GridKit::Testing::TestingResults              results;
  GridKit::Testing::SignalDelayIntegrationTests tests;
  results += tests.feedbackAndConvergence();
  results += tests.prescribedSource();
  results += tests.smoothPulseAndLongDelay();
  results += tests.rollbackAndContinuation();
  return results.summary();
}
