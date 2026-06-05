#include "GovernorHygovTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::GovernorHygovTests<double, size_t> test;

  result += test.constructionAndValidation();
  result += test.signals();
  result += test.sourceDefault();
  result += test.baseConversion();
  result += test.jsonParseAndSystemAssembly();

  return result.summary();
}
