#include "SystemModelTests.hpp"

int main()
{
  using namespace GridKit;
  using namespace GridKit::Testing;

  TestingResults   result;
  SystemModelTests test;

  result += test.signalError();

  return result.summary();
}
