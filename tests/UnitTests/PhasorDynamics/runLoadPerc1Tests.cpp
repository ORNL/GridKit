#include "LoadPerc1Tests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;
  GridKit::Testing::LoadPerc1Tests test;

  result += test.validation();
  result += test.initialization();
  result += test.initializationDomainAndFallback();
  result += test.electricalResponse();
  result += test.historyRates();
  result += test.jsonSystem();
  result += test.jacobian();

  return result.summary();
}
