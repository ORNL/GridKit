#include "ConverterRegcaTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::ConverterRegcaTests<double, size_t> test;

  result += test.constructor();
  result += test.lifecycle();
  result += test.signalVerification();
  result += test.nullBusVerification();
  result += test.jsonParseAndSystemAssembly();

  return result.summary();
}
