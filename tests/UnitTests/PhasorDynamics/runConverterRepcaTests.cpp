#include "ConverterRepcaTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::ConverterRepcaTests<double, size_t> test;

  result += test.constructionAndValidation();
  result += test.signalVerification();
  result += test.initializationAndResidual();
  result += test.residualEquations();
  result += test.jsonParseAndSystemAssembly();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif

  return result.summary();
}
