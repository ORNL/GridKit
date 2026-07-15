#include "ExciterEsdc1aTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::ExciterEsdc1aTests<double, size_t> test;

  result += test.constructor();
  result += test.zeroInitialResidual();
  result += test.blockDiagramSemantics();
  result += test.parameterValidation();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobianStructureAndValues();
#endif

  return result.summary();
}
