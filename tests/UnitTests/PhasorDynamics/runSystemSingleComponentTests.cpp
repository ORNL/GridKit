#include "SystemSingleComponentTests.hpp"

int main()
{
  using namespace GridKit;
  using namespace GridKit::Testing;

  GridKit::Testing::TestingResults                             result;
  GridKit::Testing::SystemSingleComponentTests<double, size_t> test;

  result += test.branch();
  result += test.branchTrip();
  result += test.bus();
  result += test.busFault();
  result += test.ieeet1();
  result += test.load();
  result += test.genrou();
  result += test.genClassical();
  result += test.tgov1();

  return result.summary();
}
