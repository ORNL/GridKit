#include "OptimizationTests.hpp"

int main()
{
  using namespace GridKit::Testing;
  using namespace GridKit::Testing::Optimization;

  TestingResults                         result;
  OptimizationTests<double, std::size_t> test;

  result += test.exactSparseDerivatives();

  return result.summary();
}
