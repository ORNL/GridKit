#include "PartitionedTests.hpp"

int main()
{
  GridKit::Testing::TestingResults                   result;
  GridKit::Testing::PartitionedTests<double, size_t> test;
  result += test.integration();
  return result.summary();
}
