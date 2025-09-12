#include <GridKit/Model/PowerElectronics/Resistor/Resistor.hpp>

#include "CsrJacobianTests.hpp"
#include <examples/PowerElectronics/RLCircuit/RLCircuit.hpp>

int main()
{
  using namespace GridKit;
  using namespace GridKit::Testing;
  using namespace GridKit::LinearAlgebra;

  TestingResults                   result;
  CsrJacobianTests<double, size_t> test;

  result += test.testCooVsCsrJacobian<Resistor>("Resistor", 1, 1.0);
  result += test.testCooVsCsrJacobian<Resistor>("Resistor", 1, 2.0);

  auto rl_circuit_system  = rlCircuitSystem(1.0e-8, 1.0e-8, true, 1.0, 1.0, 1.0);
  result                 += test.testSystemCooVsCsrJacobian("RLCircuit", *rl_circuit_system);

  return result.summary();
}
