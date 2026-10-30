#include <cstddef>

#include "PartitionPowerElectronicsTests.hpp"
#include <tests/UnitTests/PowerElectronics/PartitionPowerElectronicsTests.hpp>

int main()
{
  GridKit::Testing::TestingResults                                                   result;
  GridKit::PowerElectronics::Testing::PartitionPowerElectronicsTests<double, size_t> test;

  result += test.residual();
  result += test.jacobian();
  result += test.parallelResidual();

  return result.summary();
}