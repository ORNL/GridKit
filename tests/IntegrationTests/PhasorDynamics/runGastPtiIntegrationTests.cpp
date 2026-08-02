#include <GridKit/Definitions.hpp>

#include "GastPtiIntegrationTests.hpp"

int main()
{
  GridKit::Testing::TestingResults                          result;
  GridKit::Testing::GastPtiIntegrationTests<double, size_t> test;

  result += test.genrou();
  result += test.gensal();

  return result.summary();
}
