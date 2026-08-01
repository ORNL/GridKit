#include "StateDataWriterTests.hpp"

int main()
{
  GridKit::Testing::TestingResults       result;
  GridKit::Testing::StateDataWriterTests test;

  result += test.writeCompleteState();
  result += test.writeFileAndRejectFailures();

  return result.summary();
}
