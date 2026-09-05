#include "AdmittanceComponentsTests.hpp"

int main()
{
  GridKit::Testing::TestingResults                            result;
  GridKit::Testing::AdmittanceComponentsTests<double, size_t> test;
  result += test.dependentSteadyState();
  result += test.loadSteadyState();
  result += test.voltageSteadyState();
  result += test.legacySteadyState();
  result += test.validation();
  result += test.idealShortInitialization();
  result += test.standaloneStorage();
  result += test.jacobians();
  return result.summary();
}
