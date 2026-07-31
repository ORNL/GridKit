#include "BranchTests.hpp"

int main()
{
  using namespace GridKit;
  using namespace GridKit::Testing;

  GridKit::Testing::TestingResults              result;
  GridKit::Testing::BranchTests<double, size_t> test;

  result += test.constructor();
  result += test.parameterValidation();
  result += test.dataConstructorDefaults();
  result += test.signalInputs();
  result += test.parameterSetters();
  result += test.residual();
  result += test.offNominalResidual();
  result += test.jacobian();
  result += test.offNominalJacobian();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.openJacobianStructure();
#endif

  return result.summary();
}
