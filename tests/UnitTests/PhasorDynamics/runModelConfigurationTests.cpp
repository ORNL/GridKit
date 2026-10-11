#include "ModelConfigurationTests.hpp"

int main()
{
  GridKit::Testing::TestingResults                          result;
  GridKit::Testing::ModelConfigurationTests<double, size_t> test;

  result += test.parameterReader();
  result += test.configurationChecks();

  return result.summary();
}
