#include <GridKit/Definitions.hpp>

#include "PDComponentIntegrationTests.hpp"

int main()
{
  GridKit::Testing::TestingResults                              result;
  GridKit::Testing::PDComponentIntegrationTests<double, size_t> test;

  result += test.genrouGastPtiInitialization();

  return result.summary();
}
