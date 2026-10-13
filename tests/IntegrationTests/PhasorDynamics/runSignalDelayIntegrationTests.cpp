#include "SignalDelayIntegrationTests.hpp"

int main()
{
  if (!AnalysisManager::Sundials::Ida<double, size_t>::supportsStepControl())
  {
    std::cout << "Signal delay integration requires the IDAS step-check extension\n";
    return 77;
  }
  GridKit::Testing::TestingResults              results;
  GridKit::Testing::SignalDelayIntegrationTests tests;
  results += tests.feedbackAndConvergence();
  results += tests.prescribedSource();
  results += tests.smoothPulseAndLongDelay();
  results += tests.rollbackAndContinuation();
  return results.summary();
}
