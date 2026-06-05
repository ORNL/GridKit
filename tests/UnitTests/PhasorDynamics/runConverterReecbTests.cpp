#include "ConverterReecbTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::ConverterReecbTests<double, size_t> test;

  result += test.constructionAndValidation();
  result += test.reecbSignalsInitializationAndResidual();
  result += test.rejectsHalfConnectedElectricalFeedback();
  result += test.reecbCommandSignalInitialization();
  result += test.reecbElectricalFeedbackUsesMvaBase();
  result += test.reecbReferenceFallbackAtAngle();
  result += test.initializesSaturatedStartConsistently();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.reecbSaturatedVoltagePiJacobianFinite();
#endif
  result += test.jsonParseAndSystemAssembly();

  return result.summary();
}
