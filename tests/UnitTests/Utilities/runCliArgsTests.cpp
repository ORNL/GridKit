#include "CliArgsTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;
  GridKit::Testing::CliArgsTests   test;

  result += test.construction();
  result += test.simpleParse();
  result += test.positionalAndRepeatable();
  result += test.positionalAndRepeatableErrors();
  result += test.parsingErrors();

  return result.summary();
}
