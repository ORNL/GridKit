#include "ComponentCloneTests.hpp"

int main()
{
  GridKit::Testing::CircuitComponentCloneTests<double, size_t> tests;

  GridKit::Testing::TestingResults result;

  result += tests.distributedGeneratorClone();
  result += tests.microgridLineClone();
  result += tests.microgridLoadClone();
  result += tests.microgridBusDQClone();

  return result.summary();
}
