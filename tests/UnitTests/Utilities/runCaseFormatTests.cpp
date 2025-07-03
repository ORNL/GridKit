#include "CaseFormatTests.hpp"

int main()
{
  GridKit::Testing::TestingResults                  result;
  GridKit::Testing::CaseFormatTests<double, size_t> test;

  result += test.simpleParse();

  return result.summary();
}
