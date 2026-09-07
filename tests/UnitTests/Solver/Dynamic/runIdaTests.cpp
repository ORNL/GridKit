#include "IdaTests.hpp"

int main()
{
  using namespace GridKit;
  using namespace GridKit::Testing;

  GridKit::Testing::TestingResults           result;
  GridKit::Testing::IdaTests<double, size_t> test;

  result += test.acceptedHistory();
  result += test.historyStepLimits();
  result += test.maximumSteps();
  result += test.invalidTimes();
  result += test.quadratureAndAdjoint();
  result += test.evaluationFailures();
  result += test.preservesInitialState();
  result += test.callback();
  result += test.dtMonitorZero();
  result += test.dtMonitorSuppressesEpsilonFinalStep();
  result += test.fixedStep();
  result += test.suppressAlgebraicErrors();
  result += test.consistentICType();

  return result.summary();
}
