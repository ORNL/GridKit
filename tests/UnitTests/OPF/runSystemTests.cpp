#include "SystemTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;
  GridKit::Testing::OPFSystemTests tests;

  result += tests.rejectsInvalidSystems();
  result += tests.handlesCompositeTopologyAndStateCorners();
  result += tests.allocatesExactStructures();
  result += tests.evaluatesExactDerivatives();
  result += tests.writesCompatibleSolutionState();

  return result.summary();
}
