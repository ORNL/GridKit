#include "PartitionTests.hpp"

int main()
{
  GridKit::Testing::TestingResults                 result;
  GridKit::Testing::PartitionTests<double, size_t> test;

  result += test.residualMatchesIntactCase("WECC240.case.json", "WECC240.partition.json");

  return result.summary();
}
