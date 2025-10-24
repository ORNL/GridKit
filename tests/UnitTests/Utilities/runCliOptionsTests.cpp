#include "CliOptionsTests.hpp"

int main()
{
  GridKit::Testing::TestingResults  result;
  GridKit::Testing::CliOptionsTests test;

  result += test.simpleParse();

  return result.summary();
}
