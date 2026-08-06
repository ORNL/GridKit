#include <GridKit/Definitions.hpp>

#include "ComponentConnectionTests.hpp"

int main()
{
  GridKit::Testing::TestingResults                           result;
  GridKit::Testing::ComponentConnectionTests<double, size_t> test;

  result += test.genrouEsdc1a();
  result += test.genrouHygov();
  result += test.regcaRepca();

  return result.summary();
}
