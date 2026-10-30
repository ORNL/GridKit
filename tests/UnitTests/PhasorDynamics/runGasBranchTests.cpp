#include "GasBranchTests.hpp"

int main()
{
  GridKit::Testing::TestingResults                 result;
  GridKit::Testing::GasBranchTests<double, size_t> test;

  result += test.constructor();
  result += test.parameterValidation();
  result += test.residual();
  result += test.initialization();
  result += test.initializationResidual();
  result += test.massConservation();
  result += test.absoluteTolerance();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
  result += test.fixedPressureBus();
#endif

  return result.summary();
}
