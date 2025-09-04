#include "CsrJacobianTests.hpp"
#include "Model/PowerElectronics/Resistor/Resistor.hpp"

int main()
{
  using namespace GridKit;
  using namespace GridKit::Testing;
  using namespace GridKit::LinearAlgebra;

  TestingResults                   result;
  CsrJacobianTests<double, size_t> test;

  result += test.testCooVsCsrJacobian<Resistor>("Resistor", 1, 1.0);
  result += test.testCooVsCsrJacobian<Resistor>("Resistor", 1, 2.0);

  return result.summary();
}
