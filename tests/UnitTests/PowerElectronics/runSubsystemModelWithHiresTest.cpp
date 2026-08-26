#include "SubsystemModelWithHiresTest.hpp"

int main()
{
  GridKit::Testing::TestingResults                               result;
  GridKit::Testing::SubsystemModelWithHiresTests<double, size_t> test;

  result += test.residual();
  result += test.jacobian();

  return result.summary();
}
