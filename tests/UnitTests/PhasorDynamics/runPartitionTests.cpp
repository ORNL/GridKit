#include "PartitionTests.hpp"

int main()
{
  GridKit::Testing::TestingResults                 result;
  GridKit::Testing::PartitionTests<double, size_t> test;

  result += test.twoPartitions();
  result += test.multiLevel();
  result += test.optionalKeys();
  result += test.signalGroupPlacement();
  result += test.invalidFiles();
  result += test.wecc240();

  return result.summary();
}
