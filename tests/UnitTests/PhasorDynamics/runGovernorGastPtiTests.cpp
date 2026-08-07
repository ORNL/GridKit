#include <GridKit/Utilities/Logger/Logger.hpp>

#include "GovernorGastPtiTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::GovernorGastPtiTests<double, size_t> test;
  using Log                     = GridKit::Utilities::Logger;
  const auto previous_verbosity = Log::verbosity();
  Log::setVerbosity(Log::Verbosity::NONE);

  result += test.validation();
  result += test.initializationAndSignals();
  result += test.initializationDomain();
  result += test.initializationExactness();
  result += test.residualEquations();
  result += test.governorControl();
  result += test.temperatureLimiting();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif

  Log::setVerbosity(previous_verbosity);
  return result.summary();
}
