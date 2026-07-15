#include "GovernorHygovTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::GovernorHygovTests<double, size_t> test;

  result += test.constructionAndValidation();
  result += test.signals();
  result += test.sourceDefault();
  result += test.zeroTimeConstants();
  result += test.baseConversion();
  result += test.absoluteTolerance();
  result += test.prefSignal();
  result += test.prefSignalBaseConversion();
  result += test.parameterValidation();
  result += test.signalValidation();
  result += test.jsonParseAndSystemAssembly();

  return result.summary();
}
