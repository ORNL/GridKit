#include "ComponentDefinitionTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::ComponentDefinitionTests test;

  result += test.definition();
  result += test.registry();

  return result.summary();
}
