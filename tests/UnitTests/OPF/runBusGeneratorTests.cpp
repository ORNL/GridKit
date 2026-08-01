#include "BusGeneratorTests.hpp"

int main()
{
  GridKit::Testing::TestingResults    result;
  GridKit::Testing::BusGeneratorTests test;

  result += test.polymorphicSizesAndBindings();
  result += test.busInitializationAndBounds();
  result += test.boundsAndStateOutput();
  result += test.generatorInitializationBoundsAndObjective();
  result += test.generatorConstraintsAndJacobian();
  result += test.generatorOfflineBehavior();
  result += test.generatorStateOutput();

  return result.summary();
}
