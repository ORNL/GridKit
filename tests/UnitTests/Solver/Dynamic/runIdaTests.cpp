#include "IdaTests.hpp"

int main()
{
  using namespace GridKit;
  using namespace GridKit::Testing;

  GridKit::Testing::TestingResults           result;
  GridKit::Testing::IdaTests<double, size_t> test;

  result += test.callback();
  result += test.dtMonitorZero();
  result += test.dtMonitorSuppressesEpsilonFinalStep();
  result += test.fixedStep();
  result += test.suppressAlgebraicErrors();
  result += test.initType();

  return result.summary();
}
