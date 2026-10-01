#include "ControllerReecbTests.hpp"

int main()
{
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
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif

  return result.summary();
}
