#include "IpoptTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;
  GridKit::Testing::IpoptTests     test;

  result += test.reference();
  result += test.solve();

  return result.summary();
}
