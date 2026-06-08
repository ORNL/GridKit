#include "ConverterRepcaTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::ConverterRepcaTests<double, size_t> test;

  result += test.constructionAndValidation();
  result += test.repcaSignalsInitializationAndResidual();
  result += test.rejectsHalfConnectedBranchSignals();
  result += test.rejectsLineDropCompensationWithoutCurrentSignals();
  result += test.convertsSystemBaseBranchPowerToComponentBase();
  result += test.mvaZeroUsesSystemBase();
  result += test.rejectsNonSteadyConnectedReferences();
  result += test.initializationUsesBranchSignals();
  result += test.omittedFrequencyPortsDefaultToZero();
  result += test.jsonParse();

  return result.summary();
}
