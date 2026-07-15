#include "ConverterReecaTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::ConverterReecaTests<double, size_t> test;

  result += test.constructionAndValidation();
  result += test.parameterValidation();
  result += test.reecaSignalsInitializationAndResidual();
  result += test.rejectsHalfConnectedElectricalFeedback();
  result += test.reecaCommandSignalInitialization();
  result += test.reecaElectricalFeedbackUsesMvaBase();
  result += test.reecaReferenceFallbackAtAngle();
  result += test.zeroTimeConstants();
  result += test.nonnegativeCurrentLimits();
  result += test.voltageDipAndVdl();
  result += test.outputAvailability();
  result += test.allocationAndAbsoluteTolerance();
  result += test.piInitializationContract();
  result += test.priorityInitialization();
  result += test.initializesSaturatedStartConsistently();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif
  result += test.jsonParseAndSystemAssembly();

  return result.summary();
}
