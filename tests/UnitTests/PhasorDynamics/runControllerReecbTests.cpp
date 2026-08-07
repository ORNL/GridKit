#include <GridKit/Utilities/Logger/Logger.hpp>

#include "ControllerReecbTests.hpp"

int main()
{
  using Log = GridKit::Utilities::Logger;

  const auto previous_verbosity = Log::verbosity();
  Log::setVerbosity(Log::Verbosity::NONE);

  GridKit::Testing::TestingResults                       result;
  GridKit::Testing::ControllerReecbTests<double, size_t> test;

  result += test.validation();
  result += test.initializationAndSignals();
  result += test.initializationDomain();
  result += test.initializationExactness();
  result += test.residualEquations();
  result += test.selectorConfigurations();
  result += test.voltVarReferenceBase();
  result += test.reactiveControl();
  result += test.activeCurrentControl();
  result += test.dependencyTracking();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif

  Log::setVerbosity(previous_verbosity);
  return result.summary();
}
