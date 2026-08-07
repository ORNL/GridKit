#include <GridKit/Definitions.hpp>

#include "PDIntegrationTests.hpp"

int main()
{
  GridKit::Testing::TestingResults                     result;
  GridKit::Testing::PDIntegrationTests<double, size_t> test;

  result += test.twoBusBasic();
  result += test.twoBusIeeet1();
  result += test.twoBusTgov1();
  result += test.threeBusBasic();
  result += test.threeBusClassical();
  result += test.renewableControlChainRecovery();

  return result.summary();
}
