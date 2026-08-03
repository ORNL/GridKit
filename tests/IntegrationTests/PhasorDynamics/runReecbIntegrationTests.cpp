#include <GridKit/Definitions.hpp>

#include "ReecbIntegrationTests.hpp"

int main()
{
  GridKit::Testing::TestingResults                        result;
  GridKit::Testing::ReecbIntegrationTests<double, size_t> test;

  result += test.regca();
  result += test.regcaReconstructedFeedback();
  result += test.regcaLoopResponse();
  result += test.regcaLoopRecovery();

  return result.summary();
}
